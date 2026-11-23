// LosSampler.cs -- the only script. Generate UE/UAV points, cast one ray per pair,
// measure local building statistics, export everything to CSV.
//
// WHAT IT DOES (in order)
//   1. UE   : K1 points on a horizontal plane, uniformly random inside radius R1
//             around the UEs object. A point that falls inside a building is
//             re-drawn until K1 valid points exist.
//   2. UAV  : a ring of radius D around the UAVs object. At each altitude in
//             uavHeights, `uavAzimuths` points are placed evenly around the ring.
//             UAVs that land inside a building are skipped.
//             K2 = uavHeights.Length x uavAzimuths.
//   3. RAYS : one ray for every (UE, UAV) pair. green = LoS, red = blocked.
//   4. STATS: casts a grid of downward rays around the UEs object and reports the
//             local building stock (coverage, mean / std / max roof height).
//   5. CSV  : writes every ray plus the statistics.
//
// HOW "INSIDE A BUILDING" IS TESTED
//   A short ray straight up. If the point is inside a building it hits the roof.
//   Unity ignores back faces by default, and a roof seen from below IS a back
//   face, so the test temporarily enables Physics.queriesHitBackfaces. Without
//   that the upward ray would silently miss every roof.
//   The ray must be longer than the tallest local building -- see
//   `insideTestUpDistance` (default 400 m, Harbin's tallest is 336 m).
//
// SETUP
//   1. Put the city model in the scene.
//   2. Layer "Buildings"; assign every building object to it; set buildingMask
//      to that layer ONLY. If a ground/terrain mesh is in the mask, rays hit it.
//   3. Drag the UEs and UAVs empty GameObjects into the two slots. UEs sits on the
//      ground (its y is the plane the users are generated on); UAVs only supplies
//      the x/z centre of the ring -- altitudes come from uavHeights, measured above
//      the UEs plane.
//   4. Right-click the component -> "Run all".
//
// Context menus are ordered; you can also run the steps one at a time to inspect
// the result of each before continuing.

using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using UnityEngine;

public class LosSampler : MonoBehaviour
{
    [Header("Anchors -- drag your empty GameObjects here")]
    public Transform uesRoot;
    public Transform uavsRoot;

    [Header("Blocking geometry -- buildings only")]
    public LayerMask buildingMask = ~0;

    [Header("UE: horizontal plane inside radius R1")]
    public int ueCount = 200;
    public float ueRadius = 1500f;
    [Tooltip("Height of the user above the UEs plane.")]
    public float ueHeight = 1.5f;

    [Header("UAV: ring of radius D around UAVs, at each altitude")]
    [Tooltip("Horizontal distance of the ring from the UAVs object. Only its x/z is used.")]
    public float uavRingRadius = 600f;
    [Tooltip("Altitudes measured from the UEs plane (i.e. above ground).")]
    public float[] uavHeights = { 30f, 60f, 100f, 150f, 250f, 400f, 600f };
    [Tooltip("How many UAVs around the ring at each altitude. K2 = heights x azimuths.")]
    public int uavAzimuths = 12;
    public float uavAzimuthOffsetDeg = 0f;

    [Header("Inside-building test (upward ray)")]
    [Tooltip("Longer than the tallest local building. Harbin's tallest is 336 m.")]
    public float insideTestUpDistance = 400f;
    public int maxAttemptsPerPoint = 500;

    [Header("Drawing")]
    public bool drawRays = true;
    public Color losColor = new Color(0.2f, 1f, 0.3f, 1f);
    public Color nlosColor = new Color(1f, 0.25f, 0.2f, 1f);
    public float markerSize = 20f;
    [Tooltip("Gizmos are redrawn every frame; cap the count on big K1 x K2.")]
    public int maxDrawnRays = 5000;

    [Header("Local building statistics around UEs")]
    public float statsRadius = 1500f;
    [Tooltip("Downward rays per side of the square grid.")]
    public int statsGrid = 41;
    [Tooltip("Rays are cast down from this height.")]
    public float statsTopY = 500f;
    [Tooltip("A hit below this height counts as ground, not a roof.")]
    public float roofThreshold = 1.5f;

    [Header("Export")]
    public string fileName = "los_export.csv";
    public string folderName = "LOS_data";

    [Header("Random")]
    public int seed = 12345;

    [Header("Status (read-only)")]
    [SerializeField] int madeUE;
    [SerializeField] int madeUAV;
    [SerializeField] int failedUE;
    [SerializeField] int failedUAV;
    [SerializeField] int rayTotal;
    [SerializeField] int rayLos;
    [SerializeField] float pLos;
    [SerializeField] float statCoverage;
    [SerializeField] float statMeanRoof;
    [SerializeField] float statStdRoof;
    [SerializeField] float statMaxRoof;

    // cached rays for gizmos
    Vector3[] _from = new Vector3[0];
    Vector3[] _to = new Vector3[0];
    bool[] _los = new bool[0];

    // ------------------------------------------------------------------ //
    //  Inside-building test
    // ------------------------------------------------------------------ //

    /// <summary>True if a ray straight up from `p` hits a roof above it.</summary>
    bool InsideBuilding(Vector3 p)
    {
        bool old = Physics.queriesHitBackfaces;
        Physics.queriesHitBackfaces = true;
        bool hit = Physics.Raycast(p, Vector3.up, insideTestUpDistance,
                                   buildingMask, QueryTriggerInteraction.Ignore);
        Physics.queriesHitBackfaces = old;
        return hit;
    }

    // ------------------------------------------------------------------ //
    //  Random helpers -- uniform, not biased towards the centre
    // ------------------------------------------------------------------ //

    static Vector2 InDisc(System.Random rng)
    {
        double r = System.Math.Sqrt(rng.NextDouble());
        double t = rng.NextDouble() * 2.0 * System.Math.PI;
        return new Vector2((float)(r * System.Math.Cos(t)),
                           (float)(r * System.Math.Sin(t)));
    }

    // ------------------------------------------------------------------ //
    //  1. UEs
    // ------------------------------------------------------------------ //

    [ContextMenu("1. Generate UEs")]
    public void GenerateUEs()
    {
        if (uesRoot == null) { Debug.LogError("LosSampler: assign uesRoot."); return; }

        ClearChildren(uesRoot);
        var rng = new System.Random(seed);
        Vector3 c = uesRoot.position;
        madeUE = 0;
        failedUE = 0;

        for (int i = 0; i < ueCount; i++)
        {
            bool placed = false;
            for (int a = 0; a < maxAttemptsPerPoint && !placed; a++)
            {
                Vector2 d = InDisc(rng) * ueRadius;
                Vector3 p = new Vector3(c.x + d.x, c.y, c.z + d.y);
                if (InsideBuilding(p + Vector3.up * ueHeight)) continue;

                var go = new GameObject("UE_" + i.ToString("D3"));
                go.transform.SetParent(uesRoot, false);
                go.transform.position = p;
                madeUE++;
                placed = true;
            }
            if (!placed) failedUE++;
        }

        Debug.Log("LosSampler: UEs " + madeUE + "/" + ueCount +
                  " placed (failed " + failedUE + "), radius " + ueRadius + " m");
    }

    // ------------------------------------------------------------------ //
    //  2. UAVs
    // ------------------------------------------------------------------ //

    [ContextMenu("2. Generate UAVs")]
    public void GenerateUAVs()
    {
        if (uavsRoot == null) { Debug.LogError("LosSampler: assign uavsRoot."); return; }
        if (uavHeights == null || uavHeights.Length == 0)
        { Debug.LogError("LosSampler: uavHeights is empty."); return; }

        ClearChildren(uavsRoot);
        Vector3 c = uavsRoot.position;
        float baseY = uesRoot != null ? uesRoot.position.y : uavsRoot.position.y;
        int az = Mathf.Max(1, uavAzimuths);
        madeUAV = 0;
        failedUAV = 0;

        for (int hi = 0; hi < uavHeights.Length; hi++)
        {
            float h = uavHeights[hi];
            for (int a = 0; a < az; a++)
            {
                float phi = (uavAzimuthOffsetDeg * Mathf.Deg2Rad)
                            + a * (Mathf.PI * 2f / az);
                Vector3 p = new Vector3(c.x + uavRingRadius * Mathf.Cos(phi),
                                        baseY + h,
                                        c.z + uavRingRadius * Mathf.Sin(phi));

                // a UAV inside a building is not a valid aerial position; skip it
                if (InsideBuilding(p)) { failedUAV++; continue; }

                var go = new GameObject("UAV_H" + h.ToString("F0") + "_A" + a.ToString("D2"));
                go.transform.SetParent(uavsRoot, false);
                go.transform.position = p;
                madeUAV++;
            }
        }

        Debug.Log("LosSampler: UAVs " + madeUAV + " placed (" + failedUAV +
                  " skipped inside buildings), " + uavHeights.Length + " altitudes x " +
                  az + " azimuths, ring radius " + uavRingRadius + " m");
    }

    [ContextMenu("Clear UEs")]
    public void ClearUEs() { ClearChildren(uesRoot); _from = _to = new Vector3[0]; _los = new bool[0]; }

    [ContextMenu("Clear UAVs")]
    public void ClearUAVs() { ClearChildren(uavsRoot); _from = _to = new Vector3[0]; _los = new bool[0]; }

    static void ClearChildren(Transform root)
    {
        if (root == null) return;
        for (int i = root.childCount - 1; i >= 0; i--)
        {
            var go = root.GetChild(i).gameObject;
            if (Application.isPlaying) Destroy(go); else DestroyImmediate(go);
        }
    }

    // ------------------------------------------------------------------ //
    //  3. Rays
    // ------------------------------------------------------------------ //

    [ContextMenu("3. Cast rays")]
    public void CastRays()
    {
        var ues = Children(uesRoot);
        var uavs = Children(uavsRoot);
        if (ues.Count == 0 || uavs.Count == 0)
        {
            Debug.LogError("LosSampler: generate UEs and UAVs first.");
            return;
        }

        var from = new List<Vector3>(ues.Count * uavs.Count);
        var to = new List<Vector3>(ues.Count * uavs.Count);
        var los = new List<bool>(ues.Count * uavs.Count);

        for (int i = 0; i < ues.Count; i++)
        {
            Vector3 p0 = ues[i].position + Vector3.up * ueHeight;
            for (int j = 0; j < uavs.Count; j++)
            {
                Vector3 p1 = uavs[j].position;
                Vector3 dir = p1 - p0;
                float len = dir.magnitude;
                if (len < 1e-3f) continue;
                bool clear = !Physics.Raycast(p0, dir / len, len, buildingMask,
                                              QueryTriggerInteraction.Ignore);
                from.Add(p0);
                to.Add(p1);
                los.Add(clear);
            }
        }

        _from = from.ToArray();
        _to = to.ToArray();
        _los = los.ToArray();

        rayTotal = _los.Length;
        rayLos = 0;
        for (int k = 0; k < _los.Length; k++) if (_los[k]) rayLos++;
        pLos = rayTotal > 0 ? (float)rayLos / rayTotal : 0f;

        Debug.Log("LosSampler: " + rayTotal + " rays, " + rayLos + " LoS, P(LoS) = "
                  + pLos.ToString("F4"));
    }

    // ------------------------------------------------------------------ //
    //  4. Local building statistics -- downward raycast grid
    // ------------------------------------------------------------------ //

    [ContextMenu("4. Measure local building statistics")]
    public void MeasureStatistics()
    {
        if (uesRoot == null) { Debug.LogError("LosSampler: assign uesRoot."); return; }

        Vector3 c = uesRoot.position;
        int n = Mathf.Max(5, statsGrid);
        float step = 2f * statsRadius / (n - 1);
        var heights = new List<float>(n * n);
        int covered = 0;

        for (int gy = 0; gy < n; gy++)
        {
            for (int gx = 0; gx < n; gx++)
            {
                float x = c.x - statsRadius + gx * step;
                float z = c.z - statsRadius + gy * step;
                RaycastHit h;
                if (Physics.Raycast(new Vector3(x, statsTopY, z), Vector3.down,
                                    out h, statsTopY + 200f, buildingMask,
                                    QueryTriggerInteraction.Ignore))
                {
                    if (h.point.y > roofThreshold)
                    {
                        heights.Add(h.point.y);
                        covered++;
                    }
                }
            }
        }

        int total = n * n;
        statCoverage = (float)covered / total;
        statMaxRoof = 0f;
        statMeanRoof = 0f;
        statStdRoof = 0f;
        if (heights.Count > 0)
        {
            double sum = 0, sum2 = 0;
            foreach (float h in heights) { sum += h; sum2 += (double)h * h; }
            double mean = sum / heights.Count;
            statMeanRoof = (float)mean;
            statStdRoof = (float)System.Math.Sqrt(
                System.Math.Max(0.0, sum2 / heights.Count - mean * mean));
            foreach (float h in heights) if (h > statMaxRoof) statMaxRoof = h;
        }

        Debug.Log(string.Format(CultureInfo.InvariantCulture,
            "LosSampler stats within {0:F0} m: coverage {1:P1}, mean roof {2:F1} m, std {3:F1} m, max {4:F1} m  ({5} roof samples of {6})",
            statsRadius, statCoverage, statMeanRoof, statStdRoof, statMaxRoof,
            heights.Count, total));
    }

    // ------------------------------------------------------------------ //
    //  5. Export
    // ------------------------------------------------------------------ //

    [ContextMenu("5. Export CSV")]
    public void ExportCSV()
    {
        if (_from.Length == 0) { Debug.LogError("LosSampler: cast rays first."); return; }

        var ues = Children(uesRoot);
        var uavs = Children(uavsRoot);
        var sb = new StringBuilder();

        sb.AppendLine("# LosSampler export");
        sb.AppendLine(string.Format(CultureInfo.InvariantCulture,
            "# UEs origin=({0:F1},{1:F1},{2:F1})  R1={3:F0} m  K1={4}",
            uesRoot.position.x, uesRoot.position.y, uesRoot.position.z, ueRadius, madeUE));
        string heights = string.Join(";", System.Array.ConvertAll(uavHeights,
            v => v.ToString("F0", CultureInfo.InvariantCulture)));
        sb.AppendLine(string.Format(CultureInfo.InvariantCulture,
            "# UAVs ring centre=({0:F1},{1:F1},{2:F1})  D={3:F0} m  heights=[{4}]  azimuths={5}  K2={6}",
            uavsRoot.position.x, uavsRoot.position.y, uavsRoot.position.z,
            uavRingRadius, heights, uavAzimuths, madeUAV));
        sb.AppendLine(string.Format(CultureInfo.InvariantCulture,
            "# local building stats within R={0:F0} m: coverage={1:F4} meanRoof={2:F2} stdRoof={3:F2} maxRoof={4:F2}",
            statsRadius, statCoverage, statMeanRoof, statStdRoof, statMaxRoof));
        sb.AppendLine(string.Format(CultureInfo.InvariantCulture,
            "# rays={0} los={1} p_los={2:F4} seed={3}",
            rayTotal, rayLos, pLos, seed));
        sb.AppendLine("ue_i,uav_i,ue_x,ue_y,ue_z,uav_x,uav_y,uav_z,dx,dy,dz,d2d,d3d,los");

        int k = 0;
        for (int i = 0; i < ues.Count; i++)
        {
            for (int j = 0; j < uavs.Count; j++)
            {
                if (k >= _from.Length) break;
                Vector3 a = _from[k], b = _to[k];
                Vector3 d = b - a;
                sb.AppendLine(string.Format(CultureInfo.InvariantCulture,
                    "{0},{1},{2:F2},{3:F2},{4:F2},{5:F2},{6:F2},{7:F2},{8:F2},{9:F2},{10:F2},{11:F2},{12:F2},{13}",
                    i, j, a.x, a.y, a.z, b.x, b.y, b.z,
                    d.x, d.y, d.z,
                    Mathf.Sqrt(d.x * d.x + d.z * d.z), d.magnitude, _los[k] ? 1 : 0));
                k++;
            }
        }

        string projectRoot = Path.GetFullPath(Path.Combine(Application.dataPath, ".."));
        string outDir = Path.Combine(projectRoot, folderName);
        try
        {
            Directory.CreateDirectory(outDir);
            string path = Path.Combine(outDir, fileName);
            File.WriteAllText(path, sb.ToString());
            Debug.Log("LosSampler: wrote " + path);
        }
        catch (System.Exception e)
        {
            Debug.LogError("LosSampler: could not write CSV -- " + e.Message);
        }
    }

    // ------------------------------------------------------------------ //

    [ContextMenu("Run all")]
    public void RunAll()
    {
        GenerateUEs();
        GenerateUAVs();
        CastRays();
        MeasureStatistics();
        ExportCSV();
    }

    static List<Transform> Children(Transform root)
    {
        var list = new List<Transform>();
        if (root == null) return list;
        foreach (Transform c in root) list.Add(c);
        return list;
    }

    void OnDrawGizmos()
    {
        if (!drawRays || _from.Length == 0) return;
        int n = Mathf.Min(_from.Length, Mathf.Max(1, maxDrawnRays));
        for (int i = 0; i < n; i++)
        {
            Gizmos.color = _los[i] ? losColor : nlosColor;
            Gizmos.DrawLine(_from[i], _to[i]);
            if (markerSize > 0f) Gizmos.DrawSphere(_from[i], markerSize * 0.3f);
        }
    }
}
