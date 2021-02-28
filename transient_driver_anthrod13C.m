%%%%%% driver for TMI transient tracer simulations.
%
% G. Jake Gebbie, ggebbie@whoi.edu, WHOI, 27 July 2011.
% updated: 16 Aug 2016 for the ACDC Summer School, Bonne Bay, NL,
% Canada.
% Added to TMI v7, 8 Sept 2016
% Based off: G. Gebbie & P. Huybers, "The mean age of ocean waters
% inferred from radiocarbon observations: sensitivity to surface
% sources and accounting for mixing histories", submitted, JPO.
% TMI = Total Matrix Intercomparison, Gebbie & Huybers, JPO, 2010.


%% Load the TMI tendency matrix, L, and the boundary matrix, B, 
%  which satisfies dc/dt = L*c + B*q, 
%  where c is a global tracer distribution and q is a boundary flux.
load L_4deg_2012

%% List the years for which the global tracer is output.

% set t0 to be the time before which there is no anthropogenic
% d13C.
t0 = 1850;  % ? (look at atmospheric d13C to decide)
tf = 2010;  % final time: when does the atmospheric d13C record
            % end?

years = t0:tf;
% a longer example
% years = [0:1:100 110:10:1000 1025:25:2000 2100:100:5000]; 
% Note: if years only has 2 times, MATLAB will choose the output
% frequency (presumably based on the number of timesteps)
yearsmid = (years(1:end-1)+years(2:end))./2;
NY = length(years)

%% Initial conditions.  
% choose 'predefined' or 'rectangle' 
%initial_type = 'rectangle'
%initial_type = 'predefined'
initial_type = 'anthrod13C'
make_initial_conditions_anthrod13C

%% Boundary conditions.
% choose 'fixed' or 'varying' 
%boundary_type = 'fixed';
%boundary_type = 'varying';
boundary_type = 'anthrod13C';
make_boundary_conditions_anthrod13C

% Pre-allocate arrays.
C = nan(NY,Nfield);
T = nan(NY,1);

% options for the ODE solver.
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Jacobian',L);

switch boundary_type
  
case{'fixed'}
  [T,C] = ode15s(@(t,x) get_tendency(t,x,L),years,c0,options);
case{'varying'}
  tau = 1./12; % monthly restoring timescale
  [T,C] = ode15s(@(t,x) get_tendency(t,x,L,B,isfc,tsfc,Csfc,tau),years,c0,options);
end  

save transient_output C T

% This is the old method where I explicitly included checkpointing
% to diagnose the state through time. If you don't get enough
% output from the above method, let me know and we can try to
% resurrect this code.
% $$$ NN = (NY-1)./2;
% $$$ for nn = 1:NN
% $$$   nn
% $$$   ind = nn*2-1:nn*2+1
% $$$   yrs = years(ind)
% $$$   Cin  = sq(C(ind(1),:));
% $$$   tic
% $$$   [Ttmp,Cout] = ode15s(@(t,x) get_tendency(t,x,L,Cb,B,isfc,tsfc,Csfc),yrs,Cin,options);
% $$$   [Ttmp,Cout] = ode15s(@(t,x) At*x,yrs,Cin,options);
% $$$   toc
% $$$   C(ind,:) = sq(Cout);
% $$$   T(ind) = sq(Ttmp);
% $$$   %save ttdsummary C T % save intermediate values if necessary
% $$$ end
% $$$ save transient_output C T

%% For efficiency, C is 2D, but that's hard to visualize. 
%  Translate it to a 4D array: time x depth x latitude x longitude
for nn = 1:NY
    nn
    Cfield(nn,:,:,:) = vector_to_field(squeeze(C(nn,:)),it,jt,kt);
end

%% Now let plot at a given depth at a given time.
depth_plot = 500; % choose a depth in meters.
time_plot = 3; % choose a time in years.
idepth = find(DEPTH==depth_plot);
itime  = find(T==time_plot);
figure
contourf(LON,LAT,squeeze(Cfield(itime,idepth,:,:)),0:.05:1)
ylabel('latitude [deg N]')
xlabel('longitude [deg E]')
colorbar

%% Or plot a meridional section at a given longitude and time.
lon_plot = -30;
time_plot = 3; % choose a time in years.
itime  = find(T==time_plot);
if lon_plot < 0
  ilon = find(LON== 360+lon_plot)
else
  ilon = find(LON== lon_plot);
end

figure
contourf(LAT,-DEPTH,squeeze(Cfield(itime,:,:,ilon)),0:.05:1)
xlabel('latitude [deg N]')
ylabel('depth [m]')
colorbar

%% compute the global mean timeseries.
load volume_4deg
d13C_mean = (C*v)./sum(v); % volume weighted.
 












