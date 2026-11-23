% make_figures.m
% Run this command from the project root, or pass the absolute script path:
%   run(fullfile('Paper', 'figure_new', 'make_figures.m'))
% Requires the checked exports of analyze_final.m and compare_los_forms.m.
% Saves 300 dpi PNGs and vector PDFs next to this script. No model refitting.

outDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(outDir));
rawDir = fullfile(projectDir, 'LOS_data');
dataDir = fullfile(rawDir, 'final');
pairedDir = fullfile(dataDir, 'form_comparison');
models = ["free", "harbin", "hefei", "pooled", "loro"];
runs = readFigureCsv(fullfile(dataDir, 'runs.csv'), ...
    ["run", "city", "district", "held", "source_csv", "R1_m", "n_ue", "n_uav", ...
    "n_rays", "n_used_rays", "n_used_bins", "meanRoof_m", "stdRoof_m", "h_eff_m", ...
    "u50_" + models, "k_" + models, "rmse_" + models, "bias_" + models]);
curves = readFigureCsv(fullfile(dataDir, 'curves.csv'), ...
    ["run", "theta_deg", "u", "n_rays", "n_los", "used", "p_measured", ...
    "p_" + models, "residual_" + models]);
laws = readFigureCsv(fullfile(dataDir, 'laws.csv'), ...
    ["law", "n_train_runs", "n_train_rays", "A0", "A1", "A2", "log10_K0", "K1", "K2", ...
    "R2_u50", "R2_log10_k", "training_runs"]);
paired = readFigureCsv(fullfile(pairedDir, 'form_comparison_runs.csv'), ...
    ["run", "n_rays", "n_bins", "rmse_logtan"]);

[training, blind] = validateInputs(runs, curves, laws, paired);
geometry = loadGeometry(rawDir, training(training.run == "Harbin_downtown2_3000", :));
copyVerifiedPair(pairedDir, outDir);
saveFigure(measuredFigure(training, curves), outDir, 'fig2_measured_curves');
saveFigure(lawFigure(training, laws(laws.law == "Pooled", :)), outDir, 'fig3_pooled_parameters');
saveFigure(transferFigure(training, curves), outDir, 'fig4_validation_transfer');
saveFigure(blindFigure(blind, curves), outDir, 'fig5_area7_blind');
saveFigure(geometryFigure(geometry), outDir, 'fig6_sampling_geometry');
fprintf('Figures saved to %s\n', outDir);


function t = readFigureCsv(path, requiredColumns)
% training_runs uses pipes within one field. Lock the CSV delimiter and header.
t = readtable(path, 'Delimiter', ',', 'ReadVariableNames', true, ...
    'VariableNamingRule', 'preserve', 'TextType', 'string');
columns = string(t.Properties.VariableNames);
missing = setdiff(requiredColumns, columns, 'stable');
assert(isempty(missing), 'Missing columns in %s: %s. Columns read: %s.', ...
    path, strjoin(missing, ', '), strjoin(columns, ', '));
end


function [training, blind] = validateInputs(runs, curves, laws, paired)
assert(height(runs) == 25 && numel(unique(runs.run)) == 25 && nnz(runs.held) == 1 ...
    && height(curves) == 900 && all(ismember(curves.run, runs.run)), ...
    'Unexpected run or elevation-bin inventory. Run analyze_final.m.');
training = runs(runs.held == 0, :);
blind = runs(runs.held ~= 0, :);
assert(blind.run == "Harbin_area7_1500" && sum(runs.n_rays) == 3023872 ...
    && nnz(training.city == "Harbin") == 12 && nnz(training.city == "Hefei") == 12, ...
    'The input no longer matches the 25-run manuscript dataset.');
assert(height(laws) == 3 && isequal(sort(laws.law), ["Harbin"; "Hefei"; "Pooled"]), ...
    'Expected Harbin, Hefei and pooled laws.');
assert(height(paired) == 24 && isequal(sort(paired.run), sort(training.run)), ...
    'The paired form results must come from the same 24 non-blind runs.');
for city = ["Harbin", "Hefei", "Pooled"]
    l = laws(laws.law == city, :);
    tr = training;
    if city ~= "Pooled"
        tr = tr(tr.city == city, :);
    end
    assert(l.n_train_runs == height(tr) && l.n_train_rays == sum(tr.n_rays) ...
        && isequal(sort(split(l.training_runs, "|")), sort(tr.run)), ...
        '%s law training membership is inconsistent.', city);
end

for i = 1:height(runs)
    r = runs(i, :);
    c = allBins(curves, r.run);
    b = c(c.used ~= 0, :);
    assert(height(c) == 36 && max(abs(c.theta_deg - (1.25:2.5:88.75).')) < 1e-10 ...
        && isequal(c.used ~= 0, c.n_rays >= 100) && sum(b.n_rays) == r.n_used_rays ...
        && height(b) == r.n_used_bins && sum(c.n_rays) == r.n_rays ...
        && max(abs(b.p_measured - b.n_los./b.n_rays)) < 1e-10, ...
        'Inconsistent 2.5-degree bins for %s.', r.run);
    if r.held
        assert(all(isnan([r.u50_free, r.k_free, r.u50_hefei, r.k_hefei, r.u50_loro, r.k_loro])) ...
            && all(isnan([c.p_free, c.p_hefei, c.p_loro]), 'all'), ...
            'Area 7 must have only Harbin and pooled predictions.');
        models = ["harbin", "pooled"];
    else
        p = paired(paired.run == r.run, :);
        assert(p.n_rays == r.n_used_rays && p.n_bins == r.n_used_bins ...
            && abs(p.rmse_logtan - r.rmse_free) < 1e-7, ...
            'Paired comparison disagrees with the log-tangent baseline for %s.', r.run);
        models = ["free", "harbin", "hefei", "pooled", "loro"];
    end
    for model = models
        u50 = r.("u50_" + model);
        k = r.("k_" + model);
        pred = 1 ./ (1 + exp(-k .* (b.u - u50)));
        if ismember(model, ["harbin", "hefei", "pooled"])
            l = laws(lower(laws.law) == model, :);
            h = log10(r.h_eff_m);
            expectedU50 = l.A0*h + l.A1*log10(1+r.stdRoof_m) + l.A2;
            expectedK = 10^(l.log10_K0 + l.K1*h + l.K2*log10(r.meanRoof_m));
            assert(abs(u50-expectedU50) < 1e-9 && abs(k-expectedK) < 1e-9, ...
                'Frozen %s coefficients disagree with predictions for %s.', model, r.run);
        end
        residual = b.p_measured - pred;
        w = b.n_rays / sum(b.n_rays);
        assert(max(abs(pred - b.("p_" + model))) < 1e-9 ...
            && max(abs(residual - b.("residual_" + model))) < 1e-9 ...
            && abs(r.("rmse_" + model) - sqrt(sum(w .* residual.^2))) < 1e-9 ...
            && abs(r.("bias_" + model) - mean(residual)) < 1e-9, ...
            'Prediction, RMSE or signed bias disagree for %s (%s).', r.run, model);
    end
end
fprintf('Verified 25 runs, retained bins, frozen laws, paired baseline and Area 7 exclusion.\n');
end


function copyVerifiedPair(src, dst)
for ext = [".png", ".pdf"]
    source = fullfile(src, "fig7_form_comparison" + ext);
    assert(isfile(source), 'Missing %s. Run compare_los_forms.m.', source);
    [ok, message] = copyfile(source, fullfile(dst, "fig1_form_comparison" + ext), 'f');
    assert(ok, '%s', message);
end
fprintf('Copied the verified paired-form figure as fig1_form_comparison.\n');
end


function b = allBins(curves, run)
b = sortrows(curves(curves.run == run, :), 'theta_deg');
end


function b = keptBins(curves, run)
b = allBins(curves, run);
b = b(b.used ~= 0, :);
end


function fig = measuredFigure(train, curves)
[fig, t] = newFigure(2, 2, [1080, 760]);
colors = [0 114 178; 213 94 0; 0 128 103; 139 65 145; 160 112 0; 70 75 85] / 255;
markers = {'o', 's', '^', 'd', 'v', 'x'};
orderH = ["suburb", "innerCity1", "innerCity2", "innerCity3", "downtown1", "downtown2"];
orderF = "Area" + string(1:6);
letters = 'abcd';
panel = 0;
for city = ["Harbin", "Hefei"]
    order = orderH;
    if city == "Hefei"
        order = orderF;
    end
    for radius = [1500, 3000]
        panel = panel + 1;
        ax = nexttile(t);
        h = gobjects(1, numel(order));
        for j = 1:numel(order)
            r = train(train.city == city & train.R1_m == radius & train.district == order(j), :);
            assert(height(r) == 1, 'Expected one %s district %s at %d m.', city, order(j), radius);
            b = allBins(curves, r.run);
            y = b.p_measured;
            y(b.used == 0) = NaN;
            h(j) = plot(ax, b.theta_deg, y, '-', 'Marker', markers{j}, 'MarkerSize', 3.5, ...
                'Color', colors(j, :), 'LineWidth', 1.1);
        end
        probabilityAxes(ax);
        title(ax, sprintf('(%s) %s, R_1 = %d m', letters(panel), city, radius));
        legend(ax, h, districtLabel(order), 'Location', 'southoutside', 'NumColumns', 3, 'Box', 'off');
    end
end
end


function fig = lawFigure(train, pooled)
[fig, t] = newFigure(1, 2, [960, 420]);
ref = {train.u50_free, log10(train.k_free)};
pred = {train.u50_pooled, log10(train.k_pooled)};
labels = {'u_{50}', 'log_{10}(k)'};
panel = {'(a) Transition midpoint', '(b) Log slope'};
r2 = [pooled.R2_u50, pooled.R2_log10_k];
blue = [0.00, 0.38, 0.65];
orange = [0.72, 0.30, 0.10];
for j = 1:2
    ax = nexttile(t);
    x = ref{j};
    y = pred{j};
    check = 1 - sum((y-x).^2) / sum((x-mean(x)).^2);
    assert(abs(check - r2(j)) < 1e-9, 'Pooled R^2 differs from laws.csv.');
    limits = [min([x;y]), max([x;y])];
    pad = 0.08 * diff(limits);
    limits = limits + [-pad, pad];
    eqLine = plot(ax, limits, limits, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.2);
    mask = train.city == "Harbin";
    h = scatter(ax, x(mask), y(mask), 29, blue, 'o', 'filled');
    f = scatter(ax, x(~mask), y(~mask), 35, orange, '^', 'filled');
    xlim(ax, limits); ylim(ax, limits); axis(ax, 'square');
    xlabel(ax, ['Free-fit ', labels{j}]);
    ylabel(ax, ['Pooled-law ', labels{j}]);
    title(ax, sprintf('%s, R^2 = %.3f', panel{j}, r2(j)));
    legend(ax, [h, f, eqLine], {'Harbin', 'Hefei', 'Equality'}, ...
        'Location', 'southeast', 'Box', 'off');
end
end


function fig = transferFigure(train, curves)
[fig, t] = newFigure(2, 2, [1050, 750]);
h = train(train.city == "Harbin", :);
f = train(train.city == "Hefei", :);
labels = {'Pooled LORO', 'Harbin to Hefei', 'Hefei to Harbin'};
colors = [0.33 0.33 0.33; 0.00 0.38 0.65; 0.72 0.30 0.10];
rmse = {train.rmse_loro, f.rmse_harbin, h.rmse_hefei};
bias = {train.bias_loro, f.bias_harbin, h.bias_hefei};
axRmse = nexttile(t);
axBias = nexttile(t);
for j = 1:3
    plotGroup(axRmse, rmse{j}, j, colors(j,:));
    plotGroup(axBias, bias{j}, j, colors(j,:));
    fprintf('%s (n=%d): RMSE median %.8f, mean %.8f, max %.8f; mean bias %+.8f\n', ...
        labels{j}, numel(rmse{j}), median(rmse{j}), mean(rmse{j}), max(rmse{j}), mean(bias{j}));
end
set(axRmse, 'XScale', 'log');
xlabel(axRmse, 'Ray-weighted RMSE (log scale)');
title(axRmse, '(a) Run-level prediction error');
xline(axBias, 0, ':', 'Color', [0.4, 0.4, 0.4], 'HandleVisibility', 'off');
xlabel(axBias, 'Mean residual (ray-traced minus predicted)');
title(axBias, '(b) Signed bin bias');
for ax = [axRmse, axBias]
    ylim(ax, [0.55, 3.45]);
    yticks(ax, 1:3);
    yticklabels(ax, labels);
    ax.YDir = 'reverse';
end
[~, worstF] = max(f.rmse_harbin);
[~, worstH] = max(h.rmse_hefei);
ax = nexttile(t);
transferCase(ax, f(worstF,:), curves, "harbin", colors(2,:), '(c)');
ax = nexttile(t);
transferCase(ax, h(worstH,:), curves, "hefei", colors(3,:), '(d)');
end


function plotGroup(ax, data, group, color)
x = sort(data);
lag = linspace(-0.15, 0.15, numel(x));
scatter(ax, x, group + lag(:), 24, color, 'o', 'filled');
med = median(x);
plot(ax, [med med], group + [-0.20 0.20], 'k-', 'LineWidth', 1.3);
end


function transferCase(ax, r, curves, model, color, letter)
b = keptBins(curves, r.run);
h = plot(ax, b.theta_deg, b.("p_" + model), '-', 'Color', color, 'LineWidth', 1.6);
p = plotMeasured(ax, b);
probabilityAxes(ax);
title(ax, sprintf('%s %s %s, R_1 = %d m, RMSE %.3f', ...
    letter, r.city, districtLabel(r.district), r.R1_m, r.("rmse_" + model)));
sourceLabel = upper(extractBefore(model, 2)) + extractAfter(model, 1) + " law";
legend(ax, [p, h], ["Ray tracing", sourceLabel], ...
    'Location', 'southeast', 'Box', 'off');
end


function fig = blindFigure(r, curves)
[fig, t] = newFigure(1, 2, [965, 440]);
b = keptBins(curves, r.run);
blue = [0.00, 0.38, 0.65];
orange = [0.72, 0.30, 0.10];
ax = nexttile(t);
h = plot(ax, b.theta_deg, b.p_harbin, '-', 'Color', blue, 'LineWidth', 1.6);
p = plot(ax, b.theta_deg, b.p_pooled, '--', 'Color', orange, 'LineWidth', 1.6);
m = plotMeasured(ax, b);
probabilityAxes(ax);
title(ax, '(a) Area 7 blind prediction');
legend(ax, [m,h,p], {'Ray tracing', 'Harbin law', 'Pooled law'}, 'Location', 'southeast', 'Box', 'off');

ax = nexttile(t);
h = plot(ax, b.theta_deg, b.residual_harbin, 'o-', 'Color', blue, 'MarkerSize', 3.5, 'LineWidth', 1.2);
p = plot(ax, b.theta_deg, b.residual_pooled, 's--', 'Color', orange, 'MarkerSize', 3.5, 'LineWidth', 1.2);
zero = yline(ax, 0, 'k-', 'HandleVisibility', 'off');
zero.LineWidth = 0.75;
for level = [-0.01, 0.01]
    yline(ax, level, ':', 'Color', [0.55, 0.55, 0.55], 'HandleVisibility', 'off');
end
xline(ax, 31, ':', '31 deg', 'Color', [0.45, 0.45, 0.45], 'HandleVisibility', 'off');
xlim(ax, [0, 90]); xticks(ax, 0:15:90);
xlabel(ax, 'Elevation angle \theta (deg)');
ylabel(ax, 'Ray-traced minus predicted probability');
title(ax, '(b) Signed residuals by retained bin');
legend(ax, [h,p], {'Harbin law', 'Pooled law'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

for model = ["harbin", "pooled"]
    e = b.("residual_" + model);
    [~, j] = max(abs(e));
    low = b.theta_deg <= 31;
    fprintf('Area 7 %s: RMSE %.8f, bias %+.8f, worst residual %+.8f at %.2f deg, low bins outside +/-0.01: %d/%d\n', ...
        model, r.("rmse_" + model), r.("bias_" + model), e(j), b.theta_deg(j), nnz(abs(e(low)) > 0.01), nnz(low));
end
end


function g = loadGeometry(rawDir, r)
assert(height(r) == 1, 'Missing the geometry reference run Harbin_downtown2_3000.');
path = fullfile(rawDir, r.source_csv);
fid = fopen(path, 'rt');
assert(fid >= 0, 'Cannot read geometry reference export: %s', path);
closeFile = onCleanup(@() fclose(fid));
header = strings(0,1);
skip = 0;
str = fgetl(fid);
while ischar(str)
    skip = skip + 1;
    if startsWith(str, 'ue_i,uav_i,')
        break;
    end
    header(end+1, 1) = string(str); %#ok<AGROW>
    str = fgetl(fid);
end
assert(ischar(str), 'Geometry export has no ray table.');
clear closeFile
txt = char(strjoin(header, newline));
ue = regexp(txt, 'UEs origin=\(([^)]+)\)', 'tokens', 'once');
uav = regexp(txt, 'UAVs ring centre=\(([^)]+)\)\s+D=([0-9.]+)', 'tokens', 'once');
roof = regexp(txt, 'stats within R=([0-9.]+)', 'tokens', 'once');
assert(~isempty(ue) && ~isempty(uav) && ~isempty(roof), 'Incomplete sampling metadata in %s.', path);
g.ueCenter = str2double(split(string(ue{1}), ','));
g.uavCenter = str2double(split(string(uav{1}), ','));
g.ringRadius = str2double(uav{2});
g.gridHalfWidth = str2double(roof{1});
g.groundRadius = r.R1_m;
g.run = r.run;
a = readmatrix(path, 'FileType', 'text', 'Delimiter', ',', 'NumHeaderLines', skip);
assert(isequal(size(a), [r.n_rays 14]), 'Geometry export size differs from runs.csv.');
[~, iu] = unique(a(:,1), 'stable');
[~, ia] = unique(a(:,2), 'stable');
g.ues = a(iu,3:5);
g.uavs = a(ia,6:8);
assert(size(g.ues, 1) == r.n_ue && size(g.uavs, 1) == r.n_uav ...
    && all(abs(g.ues(:,2)-g.ueCenter(2)-1.50) < 0.001), 'Unexpected sampling coordinates.');
end


function fig = geometryFigure(g)
[fig, t] = newFigure(1, 2, [1060, 480]);
blue = [0.00, 0.38, 0.65];
orange = [0.72, 0.30, 0.10];
ax = nexttile(t);
origin = g.ueCenter([1,3]);
ring = g.uavCenter([1,3]);
ang = linspace(0, 2*pi, 240);
gridX = origin(1) + g.gridHalfWidth*[-1 1 1 -1 -1];
gridZ = origin(2) + g.gridHalfWidth*[-1 -1 1 1 -1];
q = plot(ax, gridX, gridZ, ':', 'Color', [0.4, 0.4, 0.4]);
d = plot(ax, origin(1)+g.groundRadius*cos(ang), origin(2)+g.groundRadius*sin(ang), 'k-');
r = plot(ax, ring(1)+g.ringRadius*cos(ang), ring(2)+g.ringRadius*sin(ang), '-', 'Color', orange);
p = scatter(ax, g.ues(:,1), g.ues(:,3), 6, blue, 'filled');
az = unique(g.uavs(:,[1,3]), 'rows');
a = scatter(ax, az(:,1), az(:,2), 22, orange, '^', 'filled');
pad = 1.08*g.gridHalfWidth;
axis(ax, 'equal');
xlim(ax, origin(1)+[-pad, pad]); ylim(ax, origin(2)+[-pad, pad]);
xlabel(ax, 'x (m)'); ylabel(ax, 'z (m)');
title(ax, '(a) Plan view of one exported run');
legend(ax, [p,d,a,r,q], {'Users', 'User disc', 'Aerial azimuths', 'Aerial ring', 'Roof-grid domain'}, ...
    'Location', 'southoutside', 'NumColumns', 3, 'Box', 'off');

ax = nexttile(t);
target = g.uavs(1,[1,3]);
column = g.uavs(all(g.uavs(:,[1,3]) == target, 2),:);
base = g.ueCenter(2);
u = g.ues(1,:);
d = hypot(target(1)-u(1), target(2)-u(3));
top = max(column(:,2))-base;
gnd = plot(ax, [0, 1.06*d], [0,0], '-', 'Color', [0.4,0.4,0.4]);
gnd.HandleVisibility = 'off';
p = plot(ax, 0, u(2)-base, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
selections = [min(column(:,2)), median(column(:,2)), max(column(:,2))];
[~, j] = min(abs(column(:,2)-selections(2)));
selections(2) = column(j,2);
colors = [orange; 0.0 0.50 0.35; blue];
h = gobjects(1,3);
labels = strings(1,3);
for j = 1:3
    z = selections(j)-base;
    h(j) = plot(ax, [0,d], [u(2)-base,z], '-o', 'MarkerIndices', 2, ...
        'Color', colors(j,:), 'LineWidth', 1.45, 'MarkerFaceColor', colors(j,:));
    labels(j) = sprintf('h = %.0f m, \\theta = %.1f deg', z, atan2d(selections(j)-u(2), d));
end
xlim(ax, [-0.03,1.06]*d); ylim(ax, [-0.03,1.06]*top);
xlabel(ax, 'Horizontal separation d (m)');
ylabel(ax, 'Height above building-base plane (m)');
title(ax, '(b) Links to one aerial azimuth');
legend(ax, [p,h], ["User height 1.50 m", labels], 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off');
fprintf('Geometry example: %s, %d ground users, %d aerial nodes, %d plan-view azimuths.\n', ...
    g.run, size(g.ues, 1), size(g.uavs, 1), size(az, 1));
end


function names = districtLabel(names)
names = replace(names, ["innerCity", "downtown", "Area"], ["Inner city ", "Downtown ", "Area "]);
names(names == "suburb") = "Suburb";
end


function p = plotMeasured(ax, b)
p = plot(ax, b.theta_deg, b.p_measured, 'o', 'Color', [0.18,0.20,0.23], ...
    'MarkerFaceColor', 'white', 'MarkerSize', 3.8, 'LineStyle', 'none');
end


function probabilityAxes(ax)
xlim(ax, [0,90]); xticks(ax, 0:15:90);
ylim(ax, [0,1]); yticks(ax, 0:0.25:1);
xlabel(ax, 'Elevation angle \theta (deg)');
ylabel(ax, 'LoS probability');
end


function [fig, t] = newFigure(nRows, nCols, canvas)
fig = figure('Visible', 'off', 'Color', 'white', 'Units', 'pixels', ...
    'Position', [50 50 canvas], 'DefaultAxesNextPlot', 'add', ...
    'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 10, ...
    'DefaultAxesTickDir', 'out', 'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on');
t = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
end


function saveFigure(fig, folder, name)
exportgraphics(fig, fullfile(folder, name + ".png"), 'Resolution', 300);
exportgraphics(fig, fullfile(folder, name + ".pdf"), 'ContentType', 'vector');
close(fig);
fprintf('Exported %s.png and %s.pdf\n', name, name);
end
