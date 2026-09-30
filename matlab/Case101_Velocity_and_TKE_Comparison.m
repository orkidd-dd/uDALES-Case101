%% Case 101 - Velocity and TKE Comparison with Experimental Data
% Comparison of uDALES simulation results with experimental measurements
% for streamwise velocity, vertical velocity, and turbulent kinetic energy.


%% 1. Load Simulation

sim = udbase(101, ...
    'H:\udales_data\u-dales-output\101');

file = 'H:\udales_data\u-dales-output\101\tdump.101.nc';


%% 2. Extract Mean Streamwise Velocity Profile

ut = sim.load_stat_t('ut');

x_target = 32;
y_target = 32;

xt = sim.xm;
yt = sim.yt;
zt = sim.zt;

% Find grid point closest to canyon centre
[~, x_index] = min(abs(xt - x_target));
[~, y_index] = min(abs(yt - y_target));

% Extract velocity at canyon centre
u_block = squeeze(ut(x_index, y_index, :, :));

% Time-average vertical velocity profile
u_profile = mean(u_block, 2, 'omitnan');
u_profile = u_profile(:);


%% 3. Normalize Vertical Velocity Profile

H = 16.5;
u_ref = 0.14;

zt_vector = zt(:);

% Select the region inside the canyon
z_h = zt_vector <= H;

z_profile = zt_vector(z_h);
u_profile_h = u_profile(z_h);

% Normalize by friction velocity
u_norm1 = [0; u_profile_h / u_ref];
z_norm1 = [0; z_profile / H];


%% 4. Extract Horizontal Vertical-Velocity Profile

xm = sim.load_stat_t('xt');
zm = sim.load_stat_t('zm');
wt = sim.load_stat_t('wt');

% Average along the y-direction
u_xz = squeeze(mean(ut, 2, 'omitnan'));
w_xz = squeeze(mean(wt, 2, 'omitnan'));

% Extract final available time step
n = size(u_xz, 3);

% Select canyon region
L = (xm >= 24 & xm <= 40);

x = xm(L);
w1 = w_xz(L, :, n);

% Select mid-canyon height
z_mid = 8.25;

[~, iz] = min(abs(zm - z_mid));

w2 = w1(:, iz);

% Normalize vertical velocity
w_norm2 = w2(:) / u_ref;

% Normalize x-coordinate relative to canyon centre
x_norm2 = (x(:) - 32) / H;

% Apply wall boundary conditions
x_norm2 = [-0.5; x_norm2; 0.5];
w_norm2 = [0; w_norm2; 0];


%% 5. Load Experimental Velocity Data

exp_u_data = readmatrix('u_k0_expdata.csv');
exp_v_data = readmatrix('v_k0.csv');

exp_u_val = exp_u_data(:, 3);
exp_z_val = exp_u_data(:, 4);

exp_x_val = exp_v_data(:, 3);
exp_v_val = exp_v_data(:, 4);


%% 6. Compare Simulated and Experimental Velocity Profiles

figure('Color', 'w', ...
       'Position', [100, 100, 800, 700]);

max_val_u = max([abs(u_norm1); abs(exp_u_val)]);
max_val_v = max([abs(w_norm2); abs(exp_v_val)]);

lim_u = ceil(max_val_u * 1.1);
lim_v = ceil(max_val_v * 1.1);


% Streamwise velocity
ax1 = axes;
hold(ax1, 'on');

p1 = plot(ax1, u_norm1, z_norm1, ...
    'k-', ...
    'LineWidth', 2, ...
    'DisplayName', 'U (Sim)');

p2 = plot(ax1, exp_u_val, exp_z_val, ...
    'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 6, ...
    'DisplayName', 'U (Exp)');

ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.YAxisLocation = 'right';
ax1.XAxisLocation = 'bottom';

xlabel(ax1, 'U / u_*', ...
    'FontSize', 18, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

ylabel(ax1, 'z / H', ...
    'FontSize', 18, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

xlim(ax1, [-lim_u, lim_u]);
ylim(ax1, [0, 1]);

grid(ax1, 'on');


% Vertical velocity
ax2 = axes('Position', ax1.Position, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'left', ...
    'Color', 'none');

hold(ax2, 'on');

p3 = plot(ax2, x_norm2, w_norm2, ...
    'r--', ...
    'LineWidth', 2, ...
    'DisplayName', 'V (Sim)');

p4 = plot(ax2, exp_x_val, exp_v_val, ...
    'r^', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 6, ...
    'DisplayName', 'V (Exp)');

ax2.XColor = 'r';
ax2.YColor = 'r';

xlabel(ax2, 'x / H', ...
    'FontSize', 18, ...
    'Color', 'r', ...
    'FontWeight', 'bold');

ylabel(ax2, 'V / u_*', ...
    'FontSize', 18, ...
    'Color', 'r', ...
    'FontWeight', 'bold');

xlim(ax2, [-0.5, 0.5]);
ylim(ax2, [-lim_v, lim_v]);


legend([p1, p2, p3, p4], ...
    'Location', 'northeast');


%% 7. Load TKE Data

TKE = sim.load_stat_xyt('tketxyc');
TKE_t = sim.load_stat_t('tketc');

time = ncread(file, 'time');

% Select averaging period
time_idx = (time >= 200 & time <= 1500);


%% 8. Extract Vertical TKE Profile

z_h = zt <= H;

z_profile = zt(z_h);

% Time-average TKE
tke_profile = mean(TKE(:, time_idx), 2, 'omitnan');

tke_profile_h = tke_profile(z_h);

% Normalize by friction velocity squared
tke_norm_black = tke_profile_h / u_ref^2;
z_norm_black = z_profile / H;


%% 9. Extract Horizontal TKE Profile

L = (xt >= 24 & xt <= 40);
x = xt(L);

% Average TKE over y-direction and selected time period
TKE_xyz = TKE_t(:, :, :, time_idx);

tke_xz = squeeze(mean(TKE_xyz, 2, 'omitnan'));
tke_xz = mean(tke_xz, 3, 'omitnan');

% Select canyon height
iz_canyon = (zt >= 0 & zt <= H);

tke_x = mean(tke_xz(L, iz_canyon), 2, 'omitnan');

% Normalize
x_norm_red = (x(:) - 32) / H;
tke_norm_red = tke_x(:) / u_ref^2;


%% 10. Load Experimental TKE Data

exp_horiz_data = readmatrix('horizontal_TKE.csv');
exp_vert_data = readmatrix('Vertical_TKE.csv');

if size(exp_horiz_data, 2) >= 4
    exp_x_red = exp_horiz_data(:, 3);
    exp_tke_red = exp_horiz_data(:, 4);
else
    exp_x_red = exp_horiz_data(:, 1);
    exp_tke_red = exp_horiz_data(:, 2);
end

if size(exp_vert_data, 2) >= 4
    exp_tke_black = exp_vert_data(:, 3);
    exp_z_black = exp_vert_data(:, 4);
else
    exp_tke_black = exp_vert_data(:, 1);
    exp_z_black = exp_vert_data(:, 2);
end


%% 11. Compare Simulated and Experimental TKE

figure('Color', 'w', ...
       'Position', [100, 100, 800, 700]);

max_tke_red = max([ ...
    max(tke_norm_red(tke_norm_red < 50)), ...
    max(exp_tke_red(exp_tke_red < 50)), ...
    1.5]);

max_tke_black = max([ ...
    max(tke_norm_black(tke_norm_black < 50)), ...
    max(exp_tke_black(exp_tke_black < 50)), ...
    3.5]);


% Horizontal TKE profile
ax1 = axes;
hold(ax1, 'on');

p1 = plot(ax1, x_norm_red, tke_norm_red, ...
    'r--', ...
    'LineWidth', 2, ...
    'DisplayName', 'TKE Horiz (Sim)');

p2 = plot(ax1, exp_x_red, exp_tke_red, ...
    'r^', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 6, ...
    'DisplayName', 'TKE Horiz (Exp)');

set(ax1, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'left', ...
    'XColor', 'r', ...
    'YColor', 'r', ...
    'XLim', [-0.5, 0.5], ...
    'YLim', [0, ceil(max_tke_red * 1.2)]);

xlabel(ax1, 'x / H', ...
    'FontSize', 18, ...
    'Color', 'r', ...
    'FontWeight', 'bold');

ylabel(ax1, 'TKE / u_*^2', ...
    'FontSize', 18, ...
    'Color', 'r', ...
    'FontWeight', 'bold');

grid(ax1, 'on');


% Vertical TKE profile
ax2 = axes('Position', ax1.Position, ...
    'XAxisLocation', 'bottom', ...
    'YAxisLocation', 'right', ...
    'Color', 'none');

hold(ax2, 'on');

p3 = plot(ax2, tke_norm_black, z_norm_black, ...
    'k-', ...
    'LineWidth', 2, ...
    'DisplayName', 'TKE Vert (Sim)');

p4 = plot(ax2, exp_tke_black, exp_z_black, ...
    'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 6, ...
    'DisplayName', 'TKE Vert (Exp)');

set(ax2, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'XLim', [0, ceil(max_tke_black * 1.2)], ...
    'YLim', [0, 1]);

xlabel(ax2, 'TKE / u_*^2', ...
    'FontSize', 18, ...
    'Color', 'k', ...
    'FontWeight', 'bold');

ylabel(ax2, 'z / H', ...
    'FontSize', 18, ...
    'Color', 'k', ...
    'FontWeight', 'bold');


legend([p1, p2, p3, p4], ...
    'Location', 'northeast', ...
    'FontSize', 12);


%% 12. TKE Profiles at Different Heights

z_targets = [0 0.5 0.75 0.9];

figure;
hold on;

for k = 1:length(z_targets)

    z_target = z_targets(k) * H;

    [~, iz] = min(abs(zt - z_target));

    tke_x = tke_xz(L, iz);

    x_norm = (x(:) - 32) / H;
    tke_norm = tke_x(:) / u_ref^2;

    plot(x_norm, tke_norm, ...
        'DisplayName', ...
        sprintf('z/H = %.2f', z_targets(k)));

end

xlabel('x/H');
ylabel('TKE/u_*^2');

legend;
grid on;
