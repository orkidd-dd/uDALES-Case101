%% Case 101 - Geometry, Velocity Field, and TKE Analysis
% Analysis of the building geometry, surface properties, mean velocity
% fields, and turbulent kinetic energy for uDALES Case 101.


%% 1. Load Simulation and Inspect Facet Types

file = 'H:\udales_data\u-dales-output\101\fielddump.101.nc';

sim = udbase(101, ...
    'H:\udales_data\u-dales-output\101');

figure;
sim.plot_fac_type();


%% 2. Surface Albedo Distribution

% Extract albedo values for each surface facet
sim.factypes.al

albedo = sim.factypes.al(sim.facs.typeid);

figure;
sim.plot_fac(albedo);
colorbar;
clim([0 0.5]);

title('Surface Albedo');

xlim([0 sim.xlen]);
ylim([0 sim.ylen]);


%% 3. Load Statistical Velocity Data

zt   = sim.load_stat_xyt('zt');
uxyt = sim.load_stat_xyt('uxyt');
time = sim.load_stat_xyt('time');


%% 4. Load Time-Averaged Velocity Field

file = 'H:\udales_data\u-dales-output\101\tdump.101.nc';

sim = udbase(101, ...
    'H:\udales_data\u-dales-output\101');

sim.tfielddump;

ut = sim.load_stat_t('ut');

xm = sim.xm;
yt = sim.yt;
zt = sim.zt;
t  = sim.load_stat_t('time');


%% 5. Define Averaging Period

% Select statistical data between 200 s and 1500 s
idx = (t > 200) & (t <= 1500);

disp('Selected statistical times:')
disp(t(idx));


%% 6. Calculate Time-Averaged Velocity Field

% Average the velocity field over the selected statistical times
ut_avg = mean(ut(:,:,:,idx), 4, 'omitnan');

% Remove solid/building cells
ut_avg(sim.Su) = NaN;


%% 7. Horizontal Velocity Field at Mid-Height

% Select the horizontal plane closest to z = 8.5 m
zloc = 8.5;

k = find(abs(zt - zloc) == min(abs(zt - zloc)), 1);

u_slice = squeeze(ut_avg(:,:,k));

% Determine symmetric colour limits
Umax = max(abs(u_slice(:)), [], 'omitnan');


figure;

pcolor(xm, yt, u_slice');
shading flat;
axis equal tight;

clim([-Umax Umax]);

colorbar;

xlabel('x (m)');
ylabel('y (m)');
title(['Streamwise velocity at z = ', ...
    num2str(zt(k), '%8.1f'), ' m']);


%% 8. Horizontal Mean Velocity Field

% Average all selected time intervals
curut = mean(ut(:,:,:,idx), 4, 'omitnan');

% Remove solid/building cells
curut(sim.Su) = NaN;

u_slice = squeeze(curut(:,:,k))';

figure;

pcolor(xm, yt, u_slice);
shading flat;
axis equal tight;

Umax = max(abs(u_slice(:)), [], 'omitnan');
clim([-Umax Umax]);

cb = colorbar;
cb.Label.String = 'Streamwise velocity, $\bar{u}$ [m/s]';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;
cb.FontSize = 13;

xlabel('$x$ [m]', 'Interpreter', 'latex');
ylabel('$y$ [m]', 'Interpreter', 'latex');

title(['$\overline{u}(x,y,z=', ...
    num2str(zt(k), '%8.1f'), 'm)$'], ...
    'Interpreter', 'latex');


%% 9. Vertical Streamwise Velocity Field

% Select the vertical centreline of the domain
j = sim.jtot/2 + 1;

u_slice_xz = squeeze(curut(:,j,:))';

figure;

pcolor(xm, zt, u_slice_xz);
shading flat;
axis equal tight;

Umax = max(abs(u_slice_xz(:)), [], 'omitnan');
clim([-Umax Umax]);

cb = colorbar;
cb.Label.String = 'Streamwise velocity, $\bar{u}$ [m/s]';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;
cb.FontSize = 13;

xlabel('$x$ [m]', 'Interpreter', 'latex');
ylabel('$z$ [m]', 'Interpreter', 'latex');

title(['$\overline{u}(x,y=', ...
    num2str(yt(j), '%8.1f'), 'm,z)$'], ...
    'Interpreter', 'latex');


%% 10. Display Simulation Geometry

sim.geom.show();
sim.geom.show_outline();



%% 12. Load Turbulent Kinetic Energy

TKE = sim.load_stat_xyt('tketxyc');


%% 13. Check TKE Data Dimensions

disp('--------------------------------------------');
disp('TKE data dimensions:');
disp(['time:    ', mat2str(size(time))]);
disp(['zt:      ', mat2str(size(zt))]);
disp(['TKE:     ', mat2str(size(TKE))]);
disp('--------------------------------------------');


%% 14. Define Normalized Vertical Coordinate

H = 16;

% Convert vertical coordinate to column vector
zt_vector = zt(:);

% Normalize height by building height
z_norm = zt_vector / H;


%% 15. Define TKE Averaging Windows

intervals = [ ...
     200  800;
     200  900;
     200 1000;
     200 1100;
     200 1200;
     200 1300;
     200 1400;
     200 1500];


%% 16. Calculate and Plot TKE Profiles

figure('Color', 'w');

hold on;
box on;
grid on;

colororder(turbo(8));

% Loop through each averaging interval
for k = 1:size(intervals,1)

    % Start and end time
    t_start = intervals(k,1);
    t_end   = intervals(k,2);

    % Select available statistical time points
    time_index1 = time >= t_start & time <= t_end;

    % Number of selected points
    n_selected1 = sum(time_index1);

    fprintf('%d-%d s: %d time points selected\n', ...
        t_start, t_end, n_selected1);

    % Average TKE over selected time points
    tke1 = mean(TKE(:,time_index1), 2, 'omitnan');

    % Plot vertical TKE profile
    plot(tke1, z_norm, ...
        'LineWidth', 1.8, ...
        'DisplayName', ...
        sprintf('%d-%d s', t_start, t_end));

end


%% 17. Format TKE Profile Figure

xlabel('TKE [m$^2$/s$^2$]', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

ylabel('$z/H$', ...
    'Interpreter', 'latex', ...
    'FontSize', 12);

title('Turbulent Kinetic Energy Profiles');

yline(1.0, '--k', 'Building height', ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off', ...
    'LabelHorizontalAlignment', 'right');

legend('Location', 'northeast'); 
