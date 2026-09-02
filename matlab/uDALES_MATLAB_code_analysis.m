% Load uDALES Case 101 data
file = 'H:\udales_data\u-dales-output\101\fielddump.101.nc';

% Load Case 101
sim = udbase(101, 'H:\udales_data\u-dales-output\101');

figure;

% Plot the surface material types
sim.plot_fac_type();

% Hide the colorbar
format_surface_plot('', 'SkipColorbar', true);

% Get albedo values
sim.factypes.al

% Assign albedo to each surface
albedo = sim.factypes.al(sim.facs.typeid);

figure

% Plot surface albedo
sim.plot_fac(albedo)

colorbar

% Set the color range from 0 to 0.5
clim([0 0.5])

title('Albedo');

% Set the x- and y-axis limits
xlim([0 sim.xlen])
ylim([0 sim.ylen])
