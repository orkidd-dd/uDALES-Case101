%% Case 101: horizontally averaged boundary-layer profiles
% Standalone script: udbase is not required.
% Source NetCDF files are not modified.

clear;
clc;
close all;

%% 1. User settings
caseDir = ...
    'D:/udales_data/u-dales-output/101_HW1_isothermal295_7200s';

file = fullfile(caseDir, 'xytdump.101.nc');

H = 16;                         % Building height [m]
requestedWindow = [1000 7200];   % Requested averaging window [s]
zMaxOverH = 4;

normalizeProfiles = false;      % false: dimensional, true: normalized
uStar = 0.0586;                 % Reference velocity if normalized [m/s]

doExport = true;
outDir = fullfile(caseDir, 'figures_boundary_layer');

assert(isfile(file), 'File not found: %s', file);
assert(H > 0, 'Building height must be positive.');
assert(requestedWindow(2) > requestedWindow(1), ...
    'Invalid averaging window.');

%% 2. Coordinates and statistical output times
time = double(ncread(file, 'time'));
time = time(:);

zt = double(ncread(file, 'zt'));
zt = zt(:);

zm = double(ncread(file, 'zm'));
zm = zm(:);

assert(all(isfinite(time)) && all(diff(time) > 0), ...
    'Output times must be finite and strictly increasing.');

timeWindow = requestedWindow;
timeWindow(2) = min(requestedWindow(2), time(end));

assert(timeWindow(1) >= 0 && ...
       timeWindow(2) > timeWindow(1), ...
       'Requested window does not overlap the available statistics.');

% Each record represents statistics accumulated over an interval
% ending at time(i). Accumulation is assumed to start at t = 0.
%
% Partial endpoint intervals use their full-interval statistics
% as an approximation.
%
% Averaging interval variances/covariances does not reconstruct
% full-window variances/covariances if interval means differ.

left = [0; time(1:end-1)];

weights = max(0, ...
    min(time, timeWindow(2)) - max(left, timeWindow(1)));

assert(sum(weights) > 0, 'No statistical records selected.');

%% 3. Read profiles on their native vertical grids
U = readIntervalMean(file, 'uxyt', 'zt', weights);

uw = readIntervalMean(file, 'upwpxyt', 'zm', weights);
reynoldsStress = -uw;

uVariance = readIntervalMean(file, 'upuptxyc', 'zt', weights);
vVariance = readIntervalMean(file, 'vpvptxyc', 'zt', weights);
wVariance = readIntervalMean(file, 'wpwptxyc', 'zt', weights);

% Three-component resolved TKE, without SGS contribution.
kResolved = 0.5 * (uVariance + vVariance + wVariance);

% Two-component TKE retained for separate PIV comparisons.
k2D = 0.5 * (uVariance + wVariance);

assert(numel(U) == numel(zt), ...
    'Velocity and zt dimensions do not match.');

assert(numel(kResolved) == numel(zt), ...
    'TKE and zt dimensions do not match.');

assert(numel(reynoldsStress) == numel(zm), ...
    'Reynolds stress and zm dimensions do not match.');

%% 4. Print averaging information
fprintf('\n============================================\n');
fprintf('Source: %s\n', file);
fprintf('Requested window: %.3f to %.3f s\n', requestedWindow);
fprintf('Effective window: %.3f to %.3f s\n', timeWindow);
fprintf('Total averaging weight: %.3f s\n', sum(weights));
fprintf('Number of statistical records: %d\n', nnz(weights));
fprintf('============================================\n');

%% 5. Select dimensional or normalized profiles
heights = {zt, zm, zt};

panelTitles = {
    '(a) Mean streamwise velocity'
    '(b) Reynolds shear stress'
    '(c) Resolved TKE'
};

if normalizeProfiles

    assert(uStar > 0, 'Reference friction velocity must be positive.');

    values = {
        U / uStar
        reynoldsStress / uStar^2
        kResolved / uStar^2
    };

    xlabels = {
        '$\langle\overline{u}\rangle_{xy}/u_*$'
        '$-\langle\overline{u''w''}\rangle_{xy}/u_*^2$'
        '$\langle TKE\rangle_{xy}/u_*^2$'
    };

    suffix = 'normalized';

else

    values = {
        U
        reynoldsStress
        kResolved
    };

    xlabels = {
        '$\langle\overline{u}\rangle_{xy}\;(\mathrm{m\,s^{-1}})$'
        '$-\langle\overline{u''w''}\rangle_{xy}\;(\mathrm{m^2\,s^{-2}})$'
        '$\langle TKE\rangle_{xy}\;(\mathrm{m^2\,s^{-2}})$'
    };

    suffix = 'dimensional';

end

%% 6. Figure formatting
figureWidthCm = 21.1;
figureHeightCm = 12;

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2 2 figureWidthCm figureHeightCm]);

set(fig, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultTextFontName', 'Arial');

layout = tiledlayout(fig, 1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = gobjects(3, 1);

for j = 1:3

    ax(j) = nexttile(layout);

    plot(ax(j), values{j}, heights{j} / H, ...
        'Color', 'k', ...
        'LineWidth', 1.5);

    hold(ax(j), 'on');

    % Building height.
    yline(ax(j), 1, '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.0);

    ylim(ax(j), [0 zMaxOverH]);
    yticks(ax(j), 0:0.5:zMaxOverH);

    title(ax(j), panelTitles{j}, ...
        'Interpreter', 'none', ...
        'FontWeight', 'normal', ...
        'FontSize', 10);

    xlabel(ax(j), xlabels{j}, ...
        'Interpreter', 'latex', ...
        'FontSize', 11);

    set(ax(j), ...
        'FontName', 'Arial', ...
        'FontSize', 9, ...
        'LineWidth', 0.8, ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'GridColor', [0.5 0.5 0.5], ...
        'GridAlpha', 0.20, ...
        'XColor', 'k', ...
        'YColor', 'k');

    if j == 1

        ylabel(ax(j), '$z/H$', ...
            'Interpreter', 'latex', ...
            'FontSize', 12);

    else

        ax(j).YTickLabel = [];

    end

    % Use current-data limits so profiles are not clipped.
    inRange = heights{j} / H >= 0 & ...
              heights{j} / H <= zMaxOverH;

    displayed = values{j}(inRange);
    displayed = displayed(isfinite(displayed));

    assert(~isempty(displayed), ...
        'No valid profile values in panel %d.', j);

    lowerLimit = min([0; displayed(:)]);
    upperLimit = max([0; displayed(:)]);

    span = upperLimit - lowerLimit;

    if span == 0
        span = 1;
    end

    if j == 3

        xlim(ax(j), [0 upperLimit + 0.08 * span]);

    else

        xlim(ax(j), ...
            [lowerLimit - 0.06 * span, ...
             upperLimit + 0.08 * span]);

    end

end

linkaxes(ax, 'y');

% No overall title or footer, matching the reference figure.
% Include averaging information in the thesis figure caption.

drawnow;

%% 7. Optional export
if doExport

    if ~isfolder(outDir)
        mkdir(outDir);
    end

    basename = fullfile(outDir, ...
        ['baseline_7200_paper_' suffix]);

    exportgraphics(fig, [basename '.png'], ...
        'Resolution', 600, ...
        'BackgroundColor', 'white');

    exportgraphics(fig, [basename '.pdf'], ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    savefig(fig, [basename '.fig']);

    % Dimensional profiles are exported regardless of display mode.
    centreTable = table( ...
        zt, zt / H, U, kResolved, k2D, ...
        'VariableNames', { ...
        'z_m', ...
        'z_over_H', ...
        'U_m_per_s', ...
        'k_res_m2_per_s2', ...
        'k_2D_m2_per_s2'});

    writetable(centreTable, ...
        fullfile(outDir, ['cell_centre_profiles_' suffix '.csv']));

    stressTable = table( ...
        zm, zm / H, reynoldsStress, ...
        'VariableNames', { ...
        'z_m', ...
        'z_over_H', ...
        'minus_uw_m2_per_s2'});

    writetable(stressTable, ...
        fullfile(outDir, ['cell_edge_stress_' suffix '.csv']));

    intervalTable = table( ...
        left, time, weights, ...
        'VariableNames', { ...
        'interval_start_s', ...
        'interval_end_s', ...
        'overlap_weight_s'});

    writetable(intervalTable, ...
        fullfile(outDir, ['averaging_intervals_' suffix '.csv']));

    save([basename '_data.mat'], ...
        'U', 'reynoldsStress', 'kResolved', 'k2D', ...
        'zt', 'zm', 'H', 'uStar', ...
        'time', 'left', 'weights', ...
        'requestedWindow', 'timeWindow', ...
        'normalizeProfiles', 'file');

    fprintf('Saved outputs to:\n%s\n', outDir);

end

%% Local function
function profile = readIntervalMean(file, varName, zName, weights)

    variableInfo = ncinfo(file, varName);
    dimensionNames = {variableInfo.Dimensions.Name};

    iz = find(strcmp(dimensionNames, zName));
    it = find(strcmp(dimensionNames, 'time'));

    assert(numel(dimensionNames) == 2 && ...
           isscalar(iz) && isscalar(it), ...
           'Unexpected dimensions for %s. Inspect with ncinfo.', ...
           varName);

    A = double(ncread(file, varName));

    % Standardize the array to height x time.
    A = permute(A, [iz it]);

    assert(size(A, 2) == numel(weights), ...
        'Time dimension mismatch for %s.', varName);

    for k = 1:numel(variableInfo.Attributes)

        attributeName = variableInfo.Attributes(k).Name;

        if any(strcmp(attributeName, ...
                {'_FillValue', 'missing_value'}))

            missingValues = ...
                double(variableInfo.Attributes(k).Value);

            A(ismember(A, missingValues)) = NaN;

        end

    end

    selected = weights > 0;
    A = A(:, selected);
    selectedWeights = weights(selected).';

    valid = isfinite(A);
    A(~valid) = 0;

    numerator = sum(A .* selectedWeights, 2);

    denominator = sum( ...
        double(valid) .* selectedWeights, 2);

    profile = numerator ./ denominator;
    profile(denominator == 0) = NaN;

end