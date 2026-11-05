function plot_baseline_convergence_three_panels(outDir)
% Adapted from the existing MATLAB boundary-layer reader.
% Sources are read only. No udbase dependency and no stale field cache.
if nargin<1, outDir=fullfile(fileparts(mfilename('fullpath')), 'outputs', 'baseline_convergence_three_panels'); end
if ~isfolder(outDir), mkdir(outDir); end
caseDirs = {'D:/udales_data/u-dales-output/101_HW1_isothermal295_7200s'};
names = {'new_0K_7200s'};
H = 16; 
startTime = 1000;
ends = [3200 4200 5200 6200 7200];
colors = [0 .447 .741; .85 .325 .098; .929 .694 .125; .494 .184 .556; .466 .674 .188];
styles = {'-', '--', '-.', ':', '-'};
% Font sizes on the wide MATLAB canvas. At a 16 cm Word figure width,
% these render at approximately 9 pt (ticks/legend) and 10 pt (labels).
tickFontSize = 20;
labelFontSize = 22;
titleFontSize = 22;
legendFontSize = 19;
allMetrics = table();

for ic = 1:numel(caseDirs)
 file = fullfile(caseDirs{ic}, 'xytdump.101.nc');
 assert(isfile(file), 'Missing source: %s', file);
 time = double(ncread(file, 'time')); time = time(:);
 z = double(ncread(file, 'zt')); z = z(:);
 zm = double(ncread(file, 'zm')); zm = zm(:);
 R = -readProfile(file, 'upwpxyt', 'zm');
 assert(all(diff(time) > 0));
 left = [0; time(1:end-1)];
 U = readProfile(file, 'uxyt', 'zt'); K = readProfile(file, 'tketxyc', 'zt');
 sumK = .5 * (readProfile(file, 'upuptxyc', 'zt') + readProfile(file, 'vpvptxyc', 'zt') + readProfile(file, 'wpwptxyc', 'zt'));
 fprintf('%s: max TKE identity error %.8g\n', names{ic}, max(abs(K - sumK), [], 'all'));
 
 numIntervals = numel(ends);
 profilesU = zeros(numel(z), numIntervals); 
 profilesK = profilesU;
 profilesR = zeros(numel(zm), numIntervals);
 ledger = table(); sourceData = table();
 
 for j = 1:numIntervals
  % Weight each interval by its overlap with the window (starting from 1000 s)
  weights = max(0, min(time, ends(j)) - max(left, startTime));
  assert(sum(weights) > 0);
  profilesU(:, j) = weightedMean(U, weights);
  profilesK(:, j) = weightedMean(K, weights);
  profilesR(:, j) = weightedMean(R, weights);
  use = weights > 0;
  ledger = [ledger; table(ends(j), nnz(use), time(find(use, 1)), time(find(use, 1, 'last')), sum(weights), ...
    'VariableNames', {'RequestedEnd_s', 'Records', 'FirstRecord_s', 'LastRecord_s', 'RepresentedDuration_s'})]; %#ok<AGROW>
  sourceData = [sourceData; table(repmat(ends(j), numel(z), 1), z, z/H, profilesU(:, j), profilesK(:, j), ...
    'VariableNames', {'RequestedEnd_s', 'z_m', 'z_over_H', 'U_m_s', 'Resolved3C_TKE_m2_s2'})]; %#ok<AGROW>
 end
 
 fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 35.64 15], 'Visible', 'off');
 % Adjusted vertical position and height since legend is now inside ax1
 ax1 = axes(fig, 'Position', [.095*24/35.64 .18 .37*24/35.64 .68]); hold(ax1, 'on');
 ax2 = axes(fig, 'Position', [.58*24/35.64 .18 .37*24/35.64 .68]); hold(ax2, 'on');
 ax3 = axes(fig, 'Position', [1.065*24/35.64 .18 .37*24/35.64 .68]); hold(ax3, 'on');
 handles = gobjects(1, numIntervals); labels = cell(1, numIntervals);
 
 for j = 1:numIntervals
  labels{j} = sprintf('%d-%d s', startTime, ends(j));
  handles(j) = plot(ax1, profilesU(:, j), z/H, 'Color', colors(j, :), ...
    'LineStyle', styles{j}, 'LineWidth', 1.7);
  plot(ax2, profilesK(:, j), z/H, 'Color', colors(j, :), ...
    'LineStyle', styles{j}, 'LineWidth', 1.7);
  plot(ax3, profilesR(:, j), zm/H, 'Color', colors(j, :), ...
    'LineStyle', styles{j}, 'LineWidth', 1.7);
 end
 
 for ax = [ax1 ax2 ax3]
     yline(ax, 1, '--', 'Color', [.35 .35 .35], ...
         'HandleVisibility', 'off');
     ylim(ax, [0 4]);
     yticks(ax, 0:0.5:4);
     grid(ax, 'on');
     box(ax, 'on');

     set(ax, 'FontName', 'Times New Roman', ...
         'FontSize', tickFontSize, 'LineWidth', 1, 'GridAlpha', .15);
     ylabel(ax, '$z/H$', 'Interpreter', 'latex', ...
         'FontSize', labelFontSize);
 end

 title(ax1, '(a) Mean velocity', ...
     'FontSize', titleFontSize, 'FontWeight', 'normal');
 title(ax2, '(b) Resolved 3C TKE', ...
     'FontSize', titleFontSize, 'FontWeight', 'normal');
 title(ax3, '(c) Reynolds shear stress', ...
     'FontSize', titleFontSize, 'FontWeight', 'normal');

 xlabel(ax1, '$U\ (\mathrm{m/s})$', ...
     'Interpreter', 'latex', 'FontSize', labelFontSize);
 xlabel(ax2, '$\mathrm{TKE}\ (\mathrm{m^2/s^2})$', ...
     'Interpreter', 'latex', 'FontSize', labelFontSize);
 xlabel(ax3, '$-\overline{u^{\prime}w^{\prime}}\ (\mathrm{m^2/s^2})$', ...
     'Interpreter', 'latex', 'FontSize', labelFontSize);

 legend(ax1, handles, labels, 'Location', 'northwest', ...
     'Box', 'off', 'FontName', 'Times New Roman', ...
     'FontSize', legendFontSize);
 
 drawnow;
 assert(abs(ax1.Position(3) - ax2.Position(3)) < 1e-10);
 rectangles = vertcat(ax1.Position, ax2.Position, ax3.Position);
 assert(all(abs(diff(rectangles(:, 2:4), 1, 1)) < 1e-10, 'all'));
 assert(abs(diff(diff(rectangles(:, 1)))) < 1e-10);
 
 stem = fullfile(outDir, [names{ic} '_velocity_TKE_stress_convergence']);
 exportgraphics(fig, [stem '_word.png'], 'Resolution', 300);
 exportgraphics(fig, [stem '_word.pdf'], 'ContentType', 'vector');
 savefig(fig, [stem '_word.fig']);
 writetable(sourceData, [stem '_profiles.csv']);
 writetable(ledger, [stem '_windows.csv']);
 save([stem '_data.mat'], 'profilesU', 'profilesK', 'profilesR', 'z', 'zm', 'time', 'ends', 'H', 'file');
 disp(ledger);
 
 lastIdx = numel(ends);
 for j = 1:(lastIdx-1)
  du = 100 * norm(profilesU(:, j) - profilesU(:, lastIdx)) / norm(profilesU(:, lastIdx));
  dk = 100 * norm(profilesK(:, j) - profilesK(:, lastIdx)) / norm(profilesK(:, lastIdx));
  allMetrics = [allMetrics; table(string(names{ic}), ends(j), du, dk, 'VariableNames', ...
    {'Case', 'RequestedEnd_s', 'Velocity_L2_difference_percent', 'TKE_L2_difference_percent'})]; %#ok<AGROW>
 end
 close(fig);
end
writetable(allMetrics, fullfile(outDir, 'cumulative_metrics.csv'));
disp(allMetrics);
fprintf('Saved MATLAB figures and data to %s\n', outDir);
end

function A = readProfile(file, name, zName)
v = ncinfo(file, name); dims = {v.Dimensions.Name};
iz = find(strcmp(dims, zName)); it = find(strcmp(dims, 'time'));
assert(numel(dims) == 2 && isscalar(iz) && isscalar(it), 'Unexpected dimensions: %s', name);
A = permute(double(ncread(file, name)), [iz it]);
for j = 1:numel(v.Attributes)
 if any(strcmp(v.Attributes(j).Name, {'_FillValue', 'missing_value'}))
  A(ismember(A, double(v.Attributes(j).Value))) = NaN;
 end
end
end

function p = weightedMean(A, w)
use = w > 0; A = A(:, use); w = w(use).'; valid = isfinite(A); A(~valid) = 0;
p = sum(A .* w, 2) ./ sum(double(valid) .* w, 2);
end


