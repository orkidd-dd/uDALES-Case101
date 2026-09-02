% Load uDALES Case 101 data from sim.plot_fac function

% Load Case 101. udbase is use to creat  a uDLASE simulation object
sim = udbase(101, 'H:\udales_data\u-dales-output\101');

figure;

% Plot the surface material types
sim.plot_fac_type();

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
xlim([0 sim.xlen]) % set the x-axis limits from 0 to the domain length in x
ylim([0 sim.ylen])
============================================================================================================
% Load Case 101 from fielddump.101.nc(Instantaneous fields)
file='H:\udales_data\u-dales-output\101\fielddump.101.nc'
sim = udbase(101,'H:\udales_data\u-dales-output\101');

% Load the streamwise velocity and Load the x, y and z coordinates use sim.load_field function
u  = sim.load_field('u');
xm = ncread(file, 'xm');
yt = ncread(file, 'yt');
zt = ncread(file, 'zt');




