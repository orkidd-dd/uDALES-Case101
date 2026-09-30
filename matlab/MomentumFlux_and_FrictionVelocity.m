%% 
% --- Configuration ---
%======================================================= Reynolds shear stress and the subgrid-scale (SGS) momentum flux.  tau_total = <u'w'> + <u''w''>=============================================================
clear
clc
file='H:\udales_data\u-dales-output\101\xytdump.101.nc'

sim = udbase(101,'H:\udales_data\u-dales-output\101');


% 2. Load statistical data (x-y averaged over time)
zt      = sim.load_stat_xyt('zt');        % Vertical grid levels (m)
upwpxyt = sim.load_stat_xyt('upwpxyt');   % Resolved flux <u'w'>
uwxyt   = sim.load_stat_xyt('uwxyt');     % Subgrid-scale (SGS) flux <u''w''>
uxyt = sim.load_stat_xyt('uxyt');
% Calculate common axis limits for consistency
z_limit = [0 max(zt)];
flux_limit = [-0.5 0];

%Dimensionless parameters

H = 16;                         % Building height (m)

u_ref = max(uxyt(:));       % Reference velocity (m/s)

zt_vector = zt(:);              % Height vector

% Normalized height
z_norm = zt_vector / H;

% Normalize flux by Uref^2
upwpxyt_norm = upwpxyt / u_ref^2;
uwxyt_norm   = uwxyt   / u_ref^2;

% Total normalized flux
total_flux_norm = upwpxyt_norm + uwxyt_norm;

figure('Name', 'Vertical Flux Analysis', 'Color', 'w');

subplot(1, 3, 1)
hold on
for k = 1 : size(upwpxyt,2)
    plot(upwpxyt(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Resolved Flux')
xlabel("u'w' / U_{ref}^2")
ylabel('z/H')
grid on

subplot(1, 3, 2)
hold on
for k = 1 : size(uwxyt,2)
    plot(uwxyt(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Subgrid Flux')
xlabel("u'w' / U_{ref}^2")
ylabel('z/H')
grid on

subplot(1, 3, 3)
hold on
total_flux = upwpxyt + uwxyt;
for k = 1 : size(total_flux,2)
    plot(total_flux(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Total Flux')
xlabel("u'w' / U_{ref}^2")
ylabel('z/H')
grid on
%%
% --- Configuration ---
clear
clc
file='H:\udales_data\u-dales-output\101\xytdump.101.nc'

sim = udbase(101,'H:\udales_data\u-dales-output\101');


% 2. Load statistical data (x-y averaged over time)
zt      = sim.load_stat_xyt('zt');        % Vertical grid levels (m)
upwpxyt = sim.load_stat_xyt('upwpxyt');   % Resolved flux <u'w'>
uwxyt   = sim.load_stat_xyt('uwxyt');     % Subgrid-scale (SGS) flux <u''w''>
uxyt = sim.load_stat_xyt('uxyt');
% Calculate common axis limits for consistency
z_limit = [0 max(zt)];
flux_limit = [-0.5 0];
%Dimensionless parameters

H = 16;                         % Building height (m)

u_fric = sqrt(-upwpxyt);   % friction velocity
u_fric_ref= max(u_fric(:))


zt_vector = zt(:);              % Height vector

% Normalized height
z_norm = zt_vector / H;

% Normalize flux by U_fric^2
upwpxyt_norm = upwpxyt / u_fric_ref^2;
uwxyt_norm   = uwxyt   / u_fric_ref^2;

% Total normalized flux
total_flux_norm = upwpxyt_norm + uwxyt_norm;

figure('Name', 'Vertical Flux Analysis', 'Color', 'w');

subplot(1, 3, 1)
hold on
for k = 1 : size(upwpxyt,2)
    plot(upwpxyt(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Resolved Flux')
xlabel("$\overline{u'w'} / u_*^2$", 'Interpreter','latex')
ylabel('z/H')
grid on

subplot(1, 3, 2)
hold on
for k = 1 : size(uwxyt,2)
    plot(uwxyt(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Subgrid Flux')
xlabel("$\overline{u''w''} / u_*^2$", 'Interpreter','latex')
ylabel('z/H')
grid on

subplot(1, 3, 3)
hold on
total_flux = upwpxyt + uwxyt;
for k = 1 : size(total_flux,2)
    plot(total_flux(:,k),zt,'LineWidth',2)
end
yline(16,'--k');
title('Total Flux')
xlabel("$(\overline{u'w'}+\overline{u''w''}) / u_*^2$", 'Interpreter','latex')
ylabel('z/H')
grid on
%%
% Case 101 - Reynolds Shear Stress Statistical Convergence
% Progressively longer averaging windows

clear;
clc;
close all;

% File path
file = 'H:\udales_data\u-dales-output\101\xytdump.101.nc';

% Load simulation
sim = udbase(101, ...
    'H:\udales_data\u-dales-output\101');

% Load statistical data
zt      = sim.load_stat_xyt('zt');
upwpxyt = sim.load_stat_xyt('upwpxyt');
uwxyt   = sim.load_stat_xyt('uwxyt');
uxyt    = sim.load_stat_xyt('uxyt');

% Load time vector
time = ncread(file,'time');

% Display dimensions
disp('--------------------------------------------');
disp('Data dimensions:');
disp(['time:    ', mat2str(size(time))]);
disp(['zt:      ', mat2str(size(zt))]);
disp(['upwpxyt: ', mat2str(size(upwpxyt))]);
disp(['uwxyt:   ', mat2str(size(uwxyt))]);
disp(['uxyt:    ', mat2str(size(uxyt))]);
disp('--------------------------------------------');

% Building height
H = 16;

% Convert vertical coordinate to column vector
zt_vector = zt(:);

% Normalized height
z_norm = zt_vector / H;

% Progressively longer averaging intervals
intervals = [
             200 800;
             200 900
             200 1000
             200 1100
             200 1200
             200 1300
             200 1400
             200 1500];

% Create figure
figure('Color','w'); 
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
    time_index = time >= t_start & time <= t_end;

    % Number of selected points
    n_selected = sum(time_index);

    fprintf('%d-%d s: %d time points selected\n', ...
        t_start, t_end, n_selected);

    % Average Reynolds stress over selected time points
    uw_profile = mean(upwpxyt(:,time_index),2,'omitnan');

    % Reynolds shear stress
    reynolds_stress = -uw_profile;

    % Plot vertical profile
    plot(reynolds_stress, z_norm, ...
        'LineWidth',1.8, ...
        'DisplayName', ...
        sprintf('%d-%d s',t_start,t_end));

end

% Building height
yline(1,'--k', ...
    'LineWidth',1.5, ...
    'DisplayName','Building height');

% Axis labels
xlabel('$-\overline{u''w''}$ [m$^2$/s$^2$]', ...
    'Interpreter','latex', ...
    'FontSize',16);

ylabel('$z/H$', ...
    'Interpreter','latex', ...
    'FontSize',16);

% Title
title('Reynolds shear stress: progressively longer averaging windows', ...
    'FontSize',15);

% Legend
legend('Location','best');

% Axis limits
ylim([0 max(z_norm)]);

% Font size
set(gca,'FontSize',13);

% Optional: save figure
% exportgraphics(gcf, ...
%     'ReynoldsStress_Convergence_Progressive.png', ...
%     'Resolution',300);
