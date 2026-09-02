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

% height for horizontal slice
zloc = 10 

% find the closest z-index to selected height
k= find (abs(zt-zloc) == min(abs(zt-zloc)),1)

% get last time step
n= size(u,4);

%extract the horizontal slice(x,y), fix z and time
u_slice = squeeze(u(:,:,k,end));

figure
% plot a 2D matrix using colors 
pcolor(xm,yt,u_slice')

shading flat
axis equal tight % Use the same scale for x and y axes.
axis tight %Adjusts the axis limits to fit the data range tightly.

% find the maximum u, exluding NaN values
Umax = max(abs(u_slice(:)), [], 'omitnan');

% Center the color scale around zero
clim([-Umax Umax])

%Shows the colors and data values.
colorbar

xlabel('x(m)')
ylabel('y(m)')
======================================================================================================
% u at one fixed time
curut = u(:,:,:,n);
% solid locations 
curut(sim.Su) = NaN;
u_slice = squeeze(curut(:,:,k))';
figure

pcolor(xm, yt, squeeze(curut(:,:,k))');
shading flat; 
axis equal tight;
Umax = max(abs(u_slice(:)), [], 'omitnan');
clim([-Umax, Umax]); %

colorbar
xlabel('$x$ [m]', 'Interpreter','latex')
ylabel('$y$ [m]', 'Interpreter','latex')
title(['u(x,y) at z = ', num2str(zt(k)), ' m'])

==========================================================================================================
%%
j = sim.jtot/2+1; % select the index in the middle of the domain span
u_slice_xz = squeeze(curut(:,j,:))';
figure

pcolor(xm, zt, squeeze(curut(:,j,:))');
shading flat; axis equal tight; colorbar
Umax = max(abs(u_slice_xz(:)), [], 'omitnan');
clim([-Umax, Umax]); % 新版本 MATLAB 可用 clim([-Umax, Umax])

colorbar
xlabel('$x$ [m]', 'Interpreter','latex')
ylabel('$z$ [m]', 'Interpreter','latex')
title(['$\overline u(x, y=', num2str(yt(j), '%8.1f'), 'm, z)$'], 'Interpreter','latex')%%
%%
sim.geom.show();
sim.geom.show_outline();




