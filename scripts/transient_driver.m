%%%%%% driver for TMI transient tracer simulations.
%
% G. Jake Gebbie, ggebbie@whoi.edu, WHOI, 27 July 2011.
% updated: 16 Aug 2016 for the ACDC Summer School, Bonne Bay, NL,
% Canada.
%
% Added to TMI v7, 8 Sept 2016
%
% Added to TMI v8, 30 Oct 2018, after ACDC Summer School, Finse,
% Norway
%
% Based off: G. Gebbie & P. Huybers, "The mean age of ocean waters
% inferred from radiocarbon observations: sensitivity to surface
% sources and accounting for mixing histories", submitted, JPO.
% TMI = Total Matrix Intercomparison, Gebbie & Huybers, JPO, 2010.

%% Setup the TMI data.
% use this MATLAB script or 
% download manually via instructions in README
read_TMI_from_google_drive

%% Choose the TMI version of interest.
% For steady-state water-mass matrix, there are four choices.
TMIproducts = {'GH2012_2x2deg','GH2012_4x4deg'}

% choose the case
TMIno = 2; % select 1, 2,
TMIversion = TMIproducts{TMIno}

%% Load the TMI tendency matrix, L, and the boundary matrix, B, 
%  which satisfies dc/dt = L*c + B*q, 
%  where c is a global tracer distribution and q is a boundary flux.
switch TMIversion
  case 'GH2012_2x2deg'
    load L_2deg_2012
  case 'GH2012_4x4deg'
    load L_4deg_2012
  otherwise
    disp('option not available')
end

% make sure mixed-layer variables are vectors.
if size(inmixlyr,2) > 1
    inmixlyr = inmixlyr';
end
if size(imixlyr,2) > 1
    imixlyr = imixlyr';
end

%% List the years for which the global tracer is output.
years = 0:2;
% a longer example
% years = [0:1:100 110:10:1000 1025:25:2000 2100:100:5000]; 
% Note: if years only has 2 times, MATLAB will choose the output
% frequency (presumably based on the number of timesteps)
yearsmid = (years(1:end-1)+years(2:end))./2;
NY = length(years)

%% Initial conditions.  
% choose 'predefined' or 'rectangle' 
initial_type = 'rectangle'
%initial_type = 'predefined'
make_initial_conditions

%% Boundary conditions.
% choose 'fixed' or 'varying' 
boundary_type = 'fixed';
%boundary_type = 'varying';
make_boundary_conditions

% Pre-allocate arrays.
C = nan(NY,Nfield);
T = nan(NY,1);

% options for the ODE solver.
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:Nfield,'Jacobian',L);

% run the forward model with t in units of years
switch boundary_type
  
case{'fixed'}
  [T,C] = ode15s(@(t,x) get_tendency(t,x,L),years,c0,options);
case{'varying'}
  tau = 1./12; % monthly restoring timescale
  isfc = find(kt==1);
 [T,C] = ode15s(@(t,x) get_tendency(t,x,L,B,isfc,tsfc,Csfc,tau),years,c0,options);
end  

% A simple save statement to a .mat file.
%save transient_output C T

%% For efficiency, C is 2D, but that's hard to visualize. 
%  Translate it to a 4D array: time x depth x latitude x longitude
clear Cfield
for nn = 1:NY
    nn
    Cfield(nn,:,:,:) = vector_to_field(sq(C(nn,:)),it,jt,kt);
end

%% Now let plot at a given depth at a given time.
depth_plot = 500; % choose a depth in meters.
time_plot = 2; % choose a time in years.
idepth = find(DEPTH==depth_plot);
itime  = find(T==time_plot);
figure
contourf(LON,LAT,sq(Cfield(itime,idepth,:,:)),0:.05:1)
ylabel('latitude [deg N]')
xlabel('longitude [deg E]')
colorbar

%% Or plot a meridional section at a given longitude and time.
lon_plot = -30;
time_plot = 2; % choose a time in years.
itime  = find(T==time_plot);
if lon_plot < 0
    lon_plot = lon_plot + 360;
end

% make sure your choice lines up on the TMI grid.
% Warning: this while loop is not foolproof!
while isempty(ilon) 
     ilon = find(LON==lon_plot);
     lon_plot = lon_plot - 1; % in case we need to continue searching
end
 
figure
contourf(LAT,-DEPTH,sq(Cfield(itime,:,:,ilon)),0:.05:1)
xlabel('latitude [deg N]')
ylabel('depth [m]')
colorbar













