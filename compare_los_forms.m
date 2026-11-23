% compare_los_forms.m
% Paired descriptive fit of two increasing, two-parameter logistic forms.
% First run analyze_final.m to produce LOS_data/final/runs.csv and curves.csv.
% This experiment uses the same retained bins and the same ray-count weights
% in both objectives. Area 7 is excluded. No morphology law enters the fit.

projectDir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
sourceDir = fullfile(projectDir, 'LOS_data', 'final');
outputDir = fullfile(sourceDir, 'form_comparison');
runs = readtable(fullfile(sourceDir, 'runs.csv'), 'TextType', 'string');
curves = readtable(fullfile(sourceDir, 'curves.csv'), 'TextType', 'string');
held = runs.held ~= 0;
assert(nnz(held) == 1 && runs.run(held) == "Harbin_area7_1500", ...
    'Unexpected holdout in runs.csv. Regenerate it with analyze_final.m.');
runs = runs(~held, :);
assert(height(runs) == 24 && nnz(runs.city == "Harbin") == 12 && nnz(runs.city == "Hefei") == 12, ...
    'Expected 12 Harbin and 12 Hefei runs after excluding area 7.');
assert(numel(unique(runs.run)) == height(runs), 'Run IDs must be unique.');

nRuns = height(runs);
fitColumns = {'n_bins', 'n_rays', 'theta50_deg', 'k_per_deg', 'u50', 'k_logtan', ...
    'rmse_theta', 'rmse_logtan', 'delta_rmse', 'mae_theta', 'mae_logtan', 'bias_theta', 'bias_logtan'};
runResult = zeros(nRuns, numel(fitColumns));
paired = cell(nRuns, 1);
maxBaselineDifference = 0;
for i = 1:nRuns
    c = sortrows(curves(curves.run == runs.run(i) & curves.used ~= 0, :), 'theta_deg');
    theta = c.theta_deg;
    u = log10(tand(theta));
    n = c.n_rays;
    y = c.p_measured;
    assert(height(c) == runs.n_used_bins(i) && sum(n) == runs.n_used_rays(i) ...
        && all(n >= 100) && all(theta > 0 & theta < 90) && numel(unique(theta)) == height(c) ...
        && all(abs(c.u-u) < 1e-10) && all(abs(y-c.n_los./n) < 1e-10), ...
        'Binned inputs disagree with the final analysis for %s.', runs.run(i));
    w = n / sum(n);

    [theta50, kTheta] = fitForm(theta, y, w, n);
    [u50, kLog] = fitForm(u, y, w, n);
    pTheta = logistic(theta, theta50, kTheta);
    pLog = logistic(u, u50, kLog);
    maxBaselineDifference = max(maxBaselineDifference, max(abs(pLog-c.p_free)));

    eTheta = y - pTheta;
    eLog = y - pLog;
    rmseTheta = sqrt(sum(w .* eTheta.^2));
    rmseLog = sqrt(sum(w .* eLog.^2));
    runResult(i, :) = [height(c), sum(n), theta50, kTheta, u50, kLog, ...
        rmseTheta, rmseLog, rmseTheta-rmseLog, sum(w.*abs(eTheta)), sum(w.*abs(eLog)), ...
        mean(eTheta), mean(eLog)];
    paired{i} = table(c.run, theta, n, c.n_los, y, pTheta, pLog, eTheta, eLog, ...
        'VariableNames', {'run', 'theta_deg', 'n_rays', 'n_los', 'p_measured', ...
        'p_theta', 'p_logtan', 'residual_theta', 'residual_logtan'});
end
assert(isfinite(maxBaselineDifference) && maxBaselineDifference < 1e-6, ...
    'Log-tangent refits differ from analyze_final.m by %.3g. Check the source analysis.', maxBaselineDifference);
summary = [runs(:, {'run', 'city', 'district', 'R1_m'}), array2table(runResult, 'VariableNames', fitColumns)];
pairedCurves = vertcat(paired{:});

fprintf('Paired two-parameter descriptive fits: theta vs log10(tan(theta))\n');
fprintf('24 runs: 12 Harbin and 12 Hefei. Area 7 is excluded.\n');
fprintf('Same 2.5 deg retained bins (>=100 rays) and same ray weighting in both fits.\n');
fprintf('RMSE and MAE are ray-weighted within run. Bias = mean(measured - predicted) over retained bins.\n');
fprintf('Delta RMSE = RMSE(theta) - RMSE(log-tangent); positive favors log-tangent.\n');
fprintf('Maximum log-tangent probability change vs analyze_final.m: %.3g\n', maxBaselineDifference);
fprintf('\n%-25s %7s %6s %12s %12s %12s\n', 'Run', 'Rays', 'Bins', 'RMSE_theta', 'RMSE_logtan', 'Delta_RMSE');
for i = 1:nRuns
    fprintf('%-25s %7d %6d %12.8f %12.8f %+12.8f\n', summary.run(i), summary.n_rays(i), ...
        summary.n_bins(i), summary.rmse_theta(i), summary.rmse_logtan(i), summary.delta_rmse(i));
end
fprintf('\nEqual-run summaries (median/mean of paired changes calculated run by run)\n');
for city = ["All", "Harbin", "Hefei"]
    if city == "All"
        g = summary;
    else
        g = summary(summary.city == city, :);
    end
    d = g.delta_rmse;
    fprintf('%s (n=%d): RMSE theta [median %.8f, mean %.8f, max %.8f]; ', ...
        city, height(g), median(g.rmse_theta), mean(g.rmse_theta), max(g.rmse_theta));
    fprintf('log-tangent [median %.8f, mean %.8f, max %.8f]\n', ...
        median(g.rmse_logtan), mean(g.rmse_logtan), max(g.rmse_logtan));
    fprintf('  Delta [median %+.8f, mean %+.8f]; lower RMSE log-tangent %d, theta %d, tie %d\n', ...
        median(d), mean(d), nnz(d > 1e-9), nnz(d < -1e-9), nnz(abs(d) <= 1e-9));
end
[~, best] = max(summary.delta_rmse);
[~, worst] = min(summary.delta_rmse);
fprintf('Largest paired gain: %s (%+.8f); smallest: %s (%+.8f).\n', ...
    summary.run(best), summary.delta_rmse(best), summary.run(worst), summary.delta_rmse(worst));

if ~isfolder(outputDir)
    mkdir(outputDir);
end
writetable(summary, fullfile(outputDir, 'form_comparison_runs.csv'));
writetable(pairedCurves, fullfile(outputDir, 'form_comparison_curves.csv'));
fig = drawComparison(summary, pairedCurves);
exportgraphics(fig, fullfile(outputDir, 'fig7_form_comparison.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outputDir, 'fig7_form_comparison.pdf'), 'ContentType', 'vector');
close(fig);
fprintf('Wrote paired CSVs and figure to %s\n', outputDir);


function p = logistic(x, midpoint, slope)
p = 1 ./ (1 + exp(-slope .* (x-midpoint)));
end


function [midpoint, slope] = fitForm(x, y, w, n)
% A weighted affine change of coordinates improves optimizer conditioning.
% It leaves the logistic family and its ray-weighted objective unchanged.
center = sum(w .* x);
scale = sqrt(sum(w .* (x-center).^2));
z = (x-center) / scale;
assert(isfinite(scale) && scale > 0 && any(y > 0 & y < 1), 'Degenerate logistic fitting input.');
empiricalLogit = log((y.*n + 0.5) ./ ((1-y).*n + 0.5));
lineFit = ([z, ones(size(z))] .* sqrt(w)) \ (empiricalLogit .* sqrt(w));
[~, half] = min(abs(y-0.5));
starts = [repmat(z(half), 4, 1), log([0.5; 1; 2; 4])];
if lineFit(1) > 0
    starts(end+1, :) = [-lineFit(2)/lineFit(1), log(lineFit(1))];
end
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-13, ...
    'MaxIter', 1000, 'MaxFunEvals', 3000);
objective = @(q) sum(w .* (y-logistic(z, q(1), exp(q(2)))).^2);
best = [NaN, NaN];
bestLoss = Inf;
for j = 1:size(starts, 1)
    [q, loss, flag] = fminsearch(objective, starts(j, :), options);
    if flag > 0 && loss < bestLoss && all(isfinite(q)) && isfinite(exp(q(2))) && exp(q(2)) > 0
        best = q;
        bestLoss = loss;
    end
end
assert(all(isfinite(best)), 'Logistic optimizer failed to converge.');
midpoint = center + scale*best(1);
slope = exp(best(2))/scale;
end


function fig = drawComparison(summary, curves)
blue = [0.00, 0.38, 0.65];
orange = [0.72, 0.30, 0.10];
fig = figure('Visible', 'off', 'Color', 'white', 'Units', 'pixels', 'Position', [50, 50, 1270, 410]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cases = ["Harbin_suburb_1500", "Harbin_downtown1_1500"];
titles = ["(a) Harbin suburb", "(b) Harbin downtown 1"];
for j = 1:2
    ax = nexttile(layout);
    hold(ax, 'on');
    q = summary(summary.run == cases(j), :);
    c = curves(curves.run == cases(j), :);
    x = linspace(min(c.theta_deg), max(c.theta_deg), 400);
    observed = plot(ax, c.theta_deg, c.p_measured, 'o', 'Color', [0.22, 0.24, 0.28], ...
        'MarkerSize', 4, 'MarkerFaceColor', 'white', 'LineStyle', 'none');
    logfit = plot(ax, x, logistic(log10(tand(x)), q.u50, q.k_logtan), '-', 'Color', blue, 'LineWidth', 1.8);
    thetafit = plot(ax, x, logistic(x, q.theta50_deg, q.k_per_deg), '--', 'Color', orange, 'LineWidth', 1.8);
    uistack(observed, 'top');
    title(ax, titles(j));
    xlabel(ax, 'Elevation angle, \theta (deg)');
    ylabel(ax, 'LoS probability');
    xlim(ax, [0, 90]); xticks(ax, 0:30:90);
    ylim(ax, [0, 1]); yticks(ax, 0:0.25:1);
    if j == 1
        legend(ax, [observed, logfit, thetafit], {'Ray tracing', 'Log-tangent fit', 'Angle fit'}, ...
            'Location', 'southeast', 'Box', 'off');
    end
end
ax = nexttile(layout);
hold(ax, 'on');
upper = ceil(max([summary.rmse_theta; summary.rmse_logtan]) / 0.02) * 0.02;
diagLine = plot(ax, [0, upper], [0, upper], ':', 'Color', [0.42, 0.44, 0.47], 'LineWidth', 1.2);
h = summary.city == "Harbin";
hh = scatter(ax, summary.rmse_theta(h), summary.rmse_logtan(h), 28, blue, 'o', 'filled');
hf = scatter(ax, summary.rmse_theta(~h), summary.rmse_logtan(~h), 34, orange, '^', 'filled');
title(ax, '(c) Paired run errors');
xlabel(ax, 'Angle-fit RMSE');
ylabel(ax, 'Log-tangent-fit RMSE');
axis(ax, 'square');
xlim(ax, [0, upper]); ylim(ax, [0, upper]);
legend(ax, [hh, hf, diagLine], {'Harbin', 'Hefei', 'Equal RMSE'}, 'Location', 'northwest', 'Box', 'off');
for ax = findall(fig, 'Type', 'axes').'
    set(ax, 'FontSize', 11, 'FontName', 'Arial', 'LineWidth', 0.9, 'TickDir', 'out', 'Box', 'off');
    grid(ax, 'on');
end
end

