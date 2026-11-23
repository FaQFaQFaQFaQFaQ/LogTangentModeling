% analyze_final.m
% Rebuild the LoS analysis exclusively from the raw LosSampler exports.
% Run from any MATLAB working directory. Uses base MATLAB.
% runs.csv has one row per run; curves.csv has one row per run and bin.
% The 'used' column identifies the bins included in fitting and scoring.

projectDir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dataDir = fullfile(projectDir, 'LOS_data');
outputDir = fullfile(dataDir, 'final');
binEdges = (0:2.5:90).';
minBinRays = 100;
expectedRays = 3040256; % Reported inventory, checked against the exports below.

[files, cities, districts, held] = listExports(dataDir);
nRuns = numel(files);
R = cell(nRuns, 1);
for i = 1:nRuns
    R{i} = readExport(fullfile(files(i).folder, files(i).name), binEdges, minBinRays);
    R{i}.held = held(i);
    R{i}.city = cities(i);
    R{i}.district = districts(i);
    R{i}.run = sprintf('%s_%s_%d', cities(i), districts(i), R{i}.R1);
    R{i}.source_csv = string(files(i).name);
end
radius = cellfun(@(r) r.R1, R);
idxH = find(cities == "Harbin" & ~held);
idxF = find(cities == "Hefei" & ~held);
idxP = find(~held);
idxA = find(held);
checkDistrictRadii(districts(idxH), radius(idxH), 'Harbin');
checkDistrictRadii(districts(idxF), radius(idxF), 'Hefei');
assert(isequal(radius(idxA), 1500), 'Expected a single Harbin area7 run at R1=1500 m.');
counts = cellfun(@(r) r.n_rays, R);

fprintf('IEEE TVT LoS analysis from raw LosSampler exports\n');
fprintf('Loaded %d runs: %d Harbin (including area7), %d Hefei; %d rays total.\n', ...
    nRuns, numel(idxH) + numel(idxA), numel(idxF), sum(counts));
fprintf('Raw rays: Harbin %d; Hefei %d; area7 %d; training pool %d.\n', ...
    sum(counts(cities == "Harbin")), sum(counts(idxF)), counts(idxA), sum(counts(idxP)));
if sum(counts) ~= expectedRays
    fprintf('Inventory discrepancy: raw %d; reported %d; raw minus reported %+d rays.\n', ...
        sum(counts), expectedRays, sum(counts) - expectedRays);
end
fprintf('R1 and morphology from headers; coverage is a fraction; h_eff = coverage*meanRoof.\n');
fprintf('theta = atan2d(dy,d2d), dy = UAV height - 1.50 m (heights from building bases).\n');
fprintf('Bins: 2.5 deg, evaluated at bin midpoints; metrics and free fits use bins with >=100 rays.\n');
fprintf('u = log10(tand(theta)); P = 1/(1+exp(-k*(u-u50))), k > 0.\n');
fprintf('Residual = measured - predicted; bias = mean(measured - predicted) over retained bins.\n');
fprintf('Positive bias means underprediction. RMSE and MAE use ray-count weights over retained bins.\n');
fprintf('Area7 is excluded from every fit, including the free fits and all LORO folds.\n');

fitParams = nan(nRuns, 2);
for i = idxP.'
    fitParams(i, :) = fitFree(R{i});
end

lh = fitLaw(R(idxH), fitParams(idxH, :), "Harbin");
lf = fitLaw(R(idxF), fitParams(idxF, :), "Hefei");
lp = fitLaw(R(idxP), fitParams(idxP, :), "Pooled");
laws = struct2table([lh; lf; lp]);

fprintf('\n=== 1. Morphology laws (ordinary least squares across training runs) ===\n');
fprintf('u50 = A0*log10(h_eff) + A1*log10(stdRoof+1) + A2\n');
fprintf('log10(k) = log10(K0) + K1*log10(h_eff) + K2*log10(meanRoof)\n');
fprintf('%-8s %3s %14s %14s %14s %14s %14s %14s %12s %12s %12s\n', ...
    'Law', 'N', 'A0', 'A1', 'A2', 'K0', 'log10_K0', 'K1', 'K2', 'R2_u50', 'R2_log10_k');
for j = 1:height(laws)
    l = laws(j, :);
    fprintf('%-8s %3d %+14.9f %+14.9f %+14.9f %14.9f %+14.9f %+14.9f %+12.9f %12.8f %12.8f\n', ...
        l.law, l.n_train_runs, l.A0, l.A1, l.A2, l.K0, l.log10_K0, l.K1, l.K2, ...
        l.R2_u50, l.R2_log10_k);
end

% Models are evaluated on the same retained bins for each target run.
% Free, Hefei, and LORO deliberately have no predictions on area7.
modelNames = ["free", "harbin", "hefei", "pooled", "loro"];
runRows = cell(nRuns, 1);
curveRows = cell(nRuns, 1);
for i = 1:nRuns
    r = R{i};
    params = nan(numel(modelNames), 2);
    params(1, :) = fitParams(i, :);
    params(2, :) = predictParams(r, lh);
    params(4, :) = predictParams(r, lp);
    if ~r.held
        params(3, :) = predictParams(r, lf);
        train = idxP(idxP ~= i);
        fold = fitLaw(R(train), fitParams(train, :), "LORO");
        params(5, :) = predictParams(r, fold);
    end
    [runRows{i}, curveRows{i}] = evaluateRun(r, params, modelNames);
end
runs = struct2table(vertcat(runRows{:}));
curves = vertcat(curveRows{:});

fprintf('\n=== 2. Per-run morphology, parameters, and residuals ===\n');
fprintf('The pooled law uses all 24 training runs; its area7 prediction is blind.\n');
fprintf('Area7 free-fit cells are NaN by protocol. Both blind predictions are printed in block 3.\n');
fprintf('%-25s %5s %7s %10s %5s %9s %8s %8s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n', ...
    'Run', 'R1_m', 'Rays', 'Used_rays', 'Bins', 'Coverage', 'MeanRoof', 'StdRoof', 'h_eff', ...
    'u50_free', 'k_free', 'RMSE_free', 'Bias_free', 'u50_pool', 'k_pool', 'RMSE_pool', 'Bias_pool');
for j = 1:height(runs)
    t = runs(j, :);
    fprintf('%-25s %5d %7d %10d %5d %9.4f %8.2f %8.2f %9.4f %+9.5f %9.5f %9.5f %+9.5f %+9.5f %9.5f %9.5f %+9.5f\n', ...
        t.run, t.R1_m, t.n_rays, t.n_used_rays, t.n_used_bins, t.coverage, t.meanRoof_m, ...
        t.stdRoof_m, t.h_eff_m, t.u50_free, t.k_free, t.rmse_free, t.bias_free, ...
        t.u50_pooled, t.k_pooled, t.rmse_pooled, t.bias_pooled);
end

fprintf('\n=== 3. Validation across runs, cities, and blind area7 ===\n');
fprintf('Run summaries give each target run equal weight. Worst bias is the signed value with largest magnitude.\n');
fprintf('Pooled LORO: 24 predictions, each from a 23-run law, leaving out the target run only.\n');
printStats(runs(idxP, :), 'loro');
fprintf('Cross-city Harbin -> Hefei: 12 source runs and 12 target runs.\n');
printStats(runs(idxF, :), 'harbin');
fprintf('Cross-city Hefei -> Harbin: 12 source runs and 12 target runs.\n');
printStats(runs(idxH, :), 'hefei');
fprintf('Blind Harbin -> area7: 12 source runs and 1 target run.\n');
printBlind(runs(idxA, :), curves, 'harbin');
fprintf('Blind Pooled -> area7: 24 source runs and 1 target run.\n');
printBlind(runs(idxA, :), curves, 'pooled');
fprintf('For each blind model, across-run median, mean, and worst equal its single-run value.\n');

if ~isfolder(outputDir)
    mkdir(outputDir);
end
writetable(runs, fullfile(outputDir, 'runs.csv'));
writetable(curves, fullfile(outputDir, 'curves.csv'));
writetable(laws, fullfile(outputDir, 'laws.csv'));
fprintf('\nWrote runs.csv, curves.csv, and laws.csv to %s\n', outputDir);


function [files, cities, districts, held] = listExports(dataDir)
harbin = dir(fullfile(dataDir, 'within_*.csv'));
hefei = dir(fullfile(dataDir, 'HeFeiData', 'HeFeiArea*.csv'));
assert(numel(harbin) == 13 && numel(hefei) == 12, ...
    'Expected 13 Harbin exports and 12 Hefei exports; found %d and %d.', numel(harbin), numel(hefei));
files = [harbin; hefei];
cities = [repmat("Harbin", numel(harbin), 1); repmat("Hefei", numel(hefei), 1)];
districts = strings(numel(files), 1);
for i = 1:numel(files)
    if cities(i) == "Harbin"
        pattern = '_(downtown[12]|innerCity[123]|suburb|area7)\.csv$';
    else
        pattern = '^HeFei(Area[1-6])_\d+\.csv$';
    end
    match = regexp(files(i).name, pattern, 'tokens', 'once');
    assert(~isempty(match), 'Unrecognized export name: %s', files(i).name);
    districts(i) = string(match{1});
end
held = cities == "Harbin" & districts == "area7";
assert(nnz(held) == 1, 'Expected exactly one Harbin area7 export.');
end


function checkDistrictRadii(districts, radii, city)
names = unique(districts);
assert(numel(names) == 6, 'Expected six training districts in %s.', city);
for name = names.'
    assert(isequal(sort(radii(districts == name)), [1500; 3000]), ...
        'Expected one 1500 m run and one 3000 m run for %s %s.', city, name);
end
end


function r = readExport(path, edges, minBinRays)
fid = fopen(path, 'rt', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot open export: %s', path);
closeFile = onCleanup(@() fclose(fid));
columnHeader = 'ue_i,uav_i,ue_x,ue_y,ue_z,uav_x,uav_y,uav_z,dx,dy,dz,d2d,d3d,los';
headerLines = {};
line = fgetl(fid);
while ischar(line) && startsWith(strtrim(line), '#')
    headerLines{end + 1} = line; %#ok<AGROW>
    line = fgetl(fid);
end
assert(ischar(line) && strcmp(strtrim(line), columnHeader), 'Unexpected CSV column header: %s', path);

header = strjoin(headerLines, newline);
r.R1 = headerNumber(header, 'R1', path);
r.coverage = headerNumber(header, 'coverage', path);
r.meanRoof = headerNumber(header, 'meanRoof', path);
r.stdRoof = headerNumber(header, 'stdRoof', path);
r.n_ue = headerNumber(header, 'K1', path);
r.n_uav = headerNumber(header, 'K2', path);
headerRays = headerNumber(header, 'rays', path);
assert(r.coverage > 0 && r.coverage <= 1 && r.meanRoof > 0 && r.stdRoof >= 0, ...
    'Invalid morphology in %s: expected fractional coverage, positive mean and nonnegative standard deviation.', path);
assert(r.n_ue == 1024 && ismember(r.n_uav, 112:120), 'Unexpected UE/UAV counts in %s.', path);
r.heff = r.coverage * r.meanRoof;

parsed = textscan(fid, repmat('%f', 1, 14), 'Delimiter', ',', ...
    'CollectOutput', true, 'ReturnOnError', false);
d = parsed{1};
r.n_rays = size(d, 1);
assert(r.n_rays == headerRays && r.n_rays == r.n_ue * r.n_uav, ...
    'Ray count mismatch in %s: read %d, header rays %d, K1*K2 %d.', ...
    path, r.n_rays, headerRays, r.n_ue * r.n_uav);
assert(all(isfinite(d), 'all'), 'Non-finite or missing ray columns in %s.', path);
ue = d(:, 1);
uav = d(:, 2);
assert(all(ue >= 0 & ue < r.n_ue & ue == fix(ue) & uav >= 0 & uav < r.n_uav & uav == fix(uav)), ...
    'UE/UAV index out of range or nonintegral in %s.', path);
pair = ue * r.n_uav + uav + 1;
assert(all(accumarray(pair, 1, [r.n_rays, 1]) == 1), 'Missing or repeated UE/UAV pairs in %s.', path);
dy = d(:, 10);
d2d = d(:, 12);
los = d(:, 14);
assert(all(d(:, 4) == 1.50 & abs(dy - (d(:, 7) - d(:, 4))) <= 0.01501), ...
    'UE height or dy disagrees with the exported vertical coordinates in %s.', path);
assert(all(dy > 0 & d2d >= 0 & (los == 0 | los == 1)), 'Invalid dy, d2d, or los values in %s.', path);

bins = discretize(atan2d(dy, d2d), edges);
nBins = numel(edges) - 1;
r.n = accumarray(bins, 1, [nBins, 1]);
r.n_los = accumarray(bins, los, [nBins, 1]);
r.p = r.n_los ./ r.n;
r.lo = edges(1:end-1);
r.hi = edges(2:end);
r.theta = (r.lo + r.hi) / 2;
r.u = log10(tand(r.theta));
r.used = r.n >= minBinRays;
assert(nnz(r.used) >= 3, 'Fewer than three retained elevation bins in %s.', path);
end


function x = headerNumber(header, key, path)
pattern = ['(?:^|\s)' key '=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?=\s|$)'];
matches = regexp(header, pattern, 'tokens');
assert(isscalar(matches), 'Expected one %s value in the header of %s.', key, path);
x = str2double(matches{1}{1});
assert(isfinite(x), 'Invalid %s value in %s.', key, path);
end


function p = logistic(u, u50, k)
p = 1 ./ (1 + exp(-k .* (u - u50)));
end


function params = fitFree(r)
assert(~r.held, 'Blind run entered free fitting: %s', r.run);
u = r.u(r.used);
p = r.p(r.used);
w = r.n(r.used) / sum(r.n(r.used));
assert(any(p > 0 & p < 1), 'No non-saturated bins in %s.', r.run);

% Smoothed empirical logits provide one initial point. The objective uses
% the original unsmoothed probabilities and ray-count weights.
n = r.n(r.used);
empiricalLogit = log((r.n_los(r.used) + 0.5) ./ (n - r.n_los(r.used) + 0.5));
lineFit = ([u, ones(size(u))] .* sqrt(w)) \ (empiricalLogit .* sqrt(w));
[~, j] = min(abs(p - 0.5));
starts = [repmat(u(j), 3, 1), log10([1; 3; 8])];
if lineFit(1) > 0
    starts(end + 1, :) = [-lineFit(2) / lineFit(1), log10(lineFit(1))];
end
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-13, ...
    'MaxIter', 1000, 'MaxFunEvals', 3000);
objective = @(q) sum(w .* (p - logistic(u, q(1), 10^q(2))).^2);
bestSSE = Inf;
params = [NaN, NaN];
for j = 1:size(starts, 1)
    [q, sse, exitflag] = fminsearch(objective, starts(j, :), options);
    k = 10^q(2);
    if exitflag > 0 && isfinite(k) && k > 0 && sse < bestSSE
        bestSSE = sse;
        params = [q(1), k];
    end
end
assert(all(isfinite(params)), 'Two-parameter optimizer failed to converge for %s.', r.run);
end


function l = fitLaw(R, fitParams, name)
assert(~any(cellfun(@(r) r.held, R)), 'Blind run entered morphology fitting.');
assert(all(isfinite(fitParams), 'all') && all(fitParams(:, 2) > 0), 'Missing free-fit parameters.');
h = cellfun(@(r) log10(r.heff), R);
s = cellfun(@(r) log10(r.stdRoof + 1), R);
hb = cellfun(@(r) log10(r.meanRoof), R);
X1 = [h, s, ones(size(h))];
X2 = [h, hb, ones(size(h))];
assert(rank(X1) == 3 && rank(X2) == 3, 'Rank-deficient morphology regression: %s.', name);
u50 = fitParams(:, 1);
logk = log10(fitParams(:, 2));
A = X1 \ u50;
B = X2 \ logk;
l = struct('law', name, 'n_train_runs', numel(R), ...
    'n_train_rays', sum(cellfun(@(r) r.n_rays, R)), ...
    'A0', A(1), 'A1', A(2), 'A2', A(3), 'K0', 10^B(3), 'log10_K0', B(3), ...
    'K1', B(1), 'K2', B(2), 'R2_u50', rSquared(u50, X1 * A), ...
    'R2_log10_k', rSquared(logk, X2 * B), ...
    'training_runs', strjoin(string(cellfun(@(r) r.run, R, 'UniformOutput', false)), '|'));
end


function value = rSquared(y, estimate)
sst = sum((y - mean(y)).^2);
value = NaN;
if sst > 0
    value = 1 - sum((y - estimate).^2) / sst;
end
end


function params = predictParams(r, l)
h = log10(r.heff);
u50 = l.A0 * h + l.A1 * log10(r.stdRoof + 1) + l.A2;
k = 10^(l.log10_K0 + l.K1 * h + l.K2 * log10(r.meanRoof));
assert(isfinite(u50) && isfinite(k) && k > 0, 'Invalid %s prediction for %s.', l.law, r.run);
params = [u50, k];
end


function [row, bins] = evaluateRun(r, params, modelNames)
row = struct('run', string(r.run), 'city', r.city, 'district', r.district, ...
    'held', r.held, 'source_csv', r.source_csv, 'R1_m', r.R1, ...
    'n_ue', r.n_ue, 'n_uav', r.n_uav, 'n_rays', r.n_rays, 'n_used_rays', sum(r.n(r.used)), ...
    'n_used_bins', nnz(r.used), 'coverage', r.coverage, 'meanRoof_m', r.meanRoof, ...
    'stdRoof_m', r.stdRoof, 'h_eff_m', r.heff);
bins = table(repmat(string(r.run), numel(r.theta), 1), r.lo, r.hi, r.theta, r.u, ...
    r.n, r.n_los, r.used, r.p, 'VariableNames', ...
    {'run', 'theta_lo_deg', 'theta_hi_deg', 'theta_deg', 'u', 'n_rays', ...
    'n_los', 'used', 'p_measured'});
w = r.n(r.used) / sum(r.n(r.used));
for j = 1:numel(modelNames)
    suffix = modelNames(j);
    u50 = params(j, 1);
    k = params(j, 2);
    pred = logistic(r.u, u50, k);
    residual = r.p - pred;
    e = residual(r.used);
    row.("u50_" + suffix) = u50;
    row.("k_" + suffix) = k;
    row.("rmse_" + suffix) = sqrt(sum(w .* e.^2));
    row.("mae_" + suffix) = sum(w .* abs(e));
    row.("bias_" + suffix) = mean(e);
    row.("max_abs_" + suffix) = max(abs(e));
    bins.("p_" + suffix) = pred;
    bins.("residual_" + suffix) = residual;
end
end


function printStats(runs, model)
fprintf('  %-12s %12s %12s %12s  %s\n', 'Metric', 'Median', 'Mean', 'Worst', 'Worst run');
for metric = ["rmse", "mae", "bias", "max_abs"]
    values = runs.(metric + "_" + model);
    [~, worst] = max(abs(values));
    fprintf('  %-12s %+12.8f %+12.8f %+12.8f  %s\n', ...
        metric, median(values), mean(values), values(worst), runs.run(worst));
end
end


function printBlind(run, curves, model)
fprintf('  %s: u50 %+.9f; k %.9f; RMSE %.8f; MAE %.8f; bias %+.8f; max |residual| %.8f\n', ...
    run.run, run.("u50_" + model), run.("k_" + model), run.("rmse_" + model), ...
    run.("mae_" + model), run.("bias_" + model), run.("max_abs_" + model));
bins = curves(curves.run == run.run & curves.used, :);
e = bins.("residual_" + model);
[~, worst] = max(abs(e));
fprintf('  Across %d retained bins: median |residual| %.8f; mean |residual| %.8f; worst %+.8f at %.2f deg.\n', ...
    height(bins), median(abs(e)), mean(abs(e)), e(worst), bins.theta_deg(worst));
end

