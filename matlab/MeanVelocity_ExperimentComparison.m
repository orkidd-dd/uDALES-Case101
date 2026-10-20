%% Data Extraction and Processing: u and w
sim = udbase(101, 'H:\udales_data\u-dales-output\101_128');

% ---------------- Parameters ----------------
time_start = 600;
time_end   = 3500;

H = 16;                   % Building height (m)
u_ref = 0.14;             % LES friction velocity (m/s)

x_target = 32;
y_target = 32;
z_mid = H/2;

% ---------------- Load data ----------------
time = sim.load_stat_t('time');
time = time(:);

ut = sim.load_stat_t('ut');
wt = sim.load_stat_t('wt');

% Coordinates corresponding to staggered velocity components
xu = sim.xm;
yu = sim.yt;
zu = sim.zt;

xw = sim.load_stat_t('xt');
zw = sim.load_stat_t('zm');

% Select statistical output times
time_mask = time >= time_start & time <= time_end;

assert(size(ut,4) == numel(time) && size(wt,4) == numel(time), ...
    'Velocity data and time dimensions do not match.');

assert(any(time_mask), ...
    'No statistical outputs found between 600 and 3500 s.');

fprintf('Selected output times: %.1f–%.1f s (%d records)\n', ...
    min(time(time_mask)), max(time(time_mask)), nnz(time_mask));

%% 1. Vertical streamwise velocity profile (black)
[~, ix] = min(abs(xu - x_target));
[~, iy] = min(abs(yu - y_target));

% Average over selected times at the canyon centre
u_profile = squeeze(mean(ut(ix, iy, :, time_mask), 4, 'omitnan'));
u_profile = u_profile(:);

zu = zu(:);
height_mask = zu > 0 & zu <= H;

u_norm1 = [0; u_profile(height_mask) / u_ref];
z_norm1 = [0; zu(height_mask) / H];

%% 2. Horizontal vertical-velocity profile (red)
% First average over selected times, then over y
w_mean = mean(wt(:,:,:,time_mask), 4, 'omitnan');
w_xz = squeeze(mean(w_mean, 2, 'omitnan'));

% Select fluid points inside the canyon
xw = xw(:);
zw = zw(:);

x_mask = xw > 24 & xw < 40;

% Interpolate to the exact mid-height
w_profile = interp1(zw, w_xz(x_mask,:).', z_mid, 'linear');
w_profile = w_profile(:);

x_norm2 = (xw(x_mask) - 32) / H;
w_norm2 = w_profile / u_ref;

% Add zero velocity at the walls
x_norm2 = [-0.5; x_norm2; 0.5];
w_norm2 = [0; w_norm2; 0];

%% 3. Experimental data
exp_u_data = readmatrix('u_k0_expdata.csv');
exp_v_data = readmatrix('v_k0.csv');

exp_u_val = exp_u_data(:,3);
exp_z_val = exp_u_data(:,4);

exp_x_val = exp_v_data(:,3);
exp_v_val = exp_v_data(:,4);

% Remove invalid rows
valid_u = isfinite(exp_u_val) & isfinite(exp_z_val);
exp_u_val = exp_u_val(valid_u);
exp_z_val = exp_z_val(valid_u);

valid_v = isfinite(exp_x_val) & isfinite(exp_v_val);
exp_x_val = exp_x_val(valid_v);
exp_v_val = exp_v_val(valid_v);

%% 4. Plot
figure('Color','w','Position',[100 100 800 700]);

max_val_u = max(abs([u_norm1; exp_u_val]), [], 'omitnan');
max_val_v = max(abs([w_norm2; exp_v_val]), [], 'omitnan');

lim_u = max(1, ceil(max_val_u * 1.1));
lim_v = max(1, ceil(max_val_v * 1.1));

% Black axes: vertical profile of streamwise velocity
ax1 = axes('Position',[0.14 0.14 0.72 0.72]);
hold(ax1,'on');

p1 = plot(ax1, u_norm1, z_norm1, 'k-', ...
    'LineWidth',2);

p2 = plot(ax1, exp_u_val, exp_z_val, 'ko', ...
    'MarkerFaceColor','k','MarkerSize',6);

ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.YAxisLocation = 'right';
ax1.XAxisLocation = 'bottom';
ax1.FontSize = 12;

xlabel(ax1,'U / u_*','FontSize',18,'FontWeight','bold');
ylabel(ax1,'z / H','FontSize',18,'FontWeight','bold');

xlim(ax1,[-lim_u lim_u]);
ylim(ax1,[0 1]);
grid(ax1,'on');

% Red axes: horizontal profile of vertical velocity
ax2 = axes('Position',ax1.Position, ...
    'XAxisLocation','top', ...
    'YAxisLocation','left', ...
    'Color','none');

hold(ax2,'on');

p3 = plot(ax2, x_norm2, w_norm2, 'r--', ...
    'LineWidth',2);

p4 = plot(ax2, exp_x_val, exp_v_val, 'r^', ...
    'MarkerFaceColor','r','MarkerSize',6);

ax2.XColor = 'r';
ax2.YColor = 'r';
ax2.FontSize = 12;

xlabel(ax2,'x / H','FontSize',18, ...
    'Color','r','FontWeight','bold');
ylabel(ax2,'V / u_*','FontSize',18, ...
    'Color','r','FontWeight','bold');

xlim(ax2,[-0.5 0.5]);
ylim(ax2,[-lim_v lim_v]);

legend(ax2,[p1 p2 p3 p4], ...
    {'U (Sim)','U (Exp)','V (Sim)','V (Exp)'}, ...
    'Location','northeast');
xlim([-10 16]);
ylim([0 1.2]);
