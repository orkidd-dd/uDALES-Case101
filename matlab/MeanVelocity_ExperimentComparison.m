%% Case 101 - Mean Velocity Profile Comparison
% Comparison of the uDALES mean streamwise velocity profile
% with digitized experimental data.


%% 1. Load Experimental Data

exp_file = 'exp_data_velocity_profiles_digitized.xlsx';

exp_raw = readmatrix(exp_file);

% Column 1: normalized height z/H
% Column 2: normalized velocity u/u*
z_exp = exp_raw(:, 1);
u_exp = exp_raw(:, 2);


%% 2. Load uDALES Simulation Data

file = 'H:\udales_data\u-dales-output\101\tdump.101.nc';

sim = udbase(101, ...
    'H:\udales_data\u-dales-output\101');

ut = sim.load_stat_t('ut');

time_tdump = ncread(file, 'time');

xt = sim.xt;
yt = sim.yt;
zt = sim.zt;


%% 3. Extract Velocity Profile at Canyon Centre

% Canyon centre coordinates
x_target = 32;
y_target = 32;

% Find the nearest grid points
[~, x_index] = min(abs(xt - x_target));
[~, y_index] = min(abs(yt - y_target));

% Extract streamwise velocity at the canyon centre
u_block = squeeze(ut(x_index, y_index, :, :));


%% 4. Define Normalization and Averaging Period

H = 16;                 % Building height [m]
u_star = 0.14;          % Friction velocity [m/s]

zt_vector = zt(:);

% Select averaging period: 200-1500 s
idx1 = time_tdump >= 200 & time_tdump <= 1500;


%% 5. Calculate Time-Averaged Velocity Profile

% Average over the selected statistical time points
u_mean = mean(u_block(:, idx1), 2, 'omitnan');

% Normalize velocity by friction velocity
u_norm = u_mean / u_star;

% Normalize height by building height
z_norm = zt_vector / H;

% Add no-slip boundary condition at the ground
u_udales = [0; u_norm(:)];
z_udales = [0; z_norm(:)];


%% 6. Compare uDALES and Experimental Profiles

figure('Color', 'w');

hold on;

% uDALES
plot(u_udales, z_udales, '-', ...
    'LineWidth', 2.0, ...
    'Color', [0.2 0.6 0.8], ...
    'DisplayName', 'uDALES');

% Experiment
plot(u_exp, z_exp, 'ro-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'r', ...
    'DisplayName', 'Experiment');


%% 7. Mark Building Height

yline(1.0, '--k', 'Building height', ...
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');


%% 8. Format Figure

grid on;
box on;

xlabel('$u / u_*$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

ylabel('$z / H$', ...
    'Interpreter', 'latex', ...
    'FontSize', 16);

title('Mean Streamwise Velocity Profile', ...
    'FontSize', 15);

legend('Location', 'best');

set(gca, 'FontSize', 13);

xlim([-10 16]);
ylim([0 1.2]);
