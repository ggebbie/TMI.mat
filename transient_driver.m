%%%%%% driver for TMI transient tracer simulations.
%
% G. Jake Gebbie, ggebbie@whoi.edu, WHOI, 27 July 2011.
% Based off: G. Gebbie & P. Huybers, "The mean age of ocean waters
% inferred from radiocarbon observations: sensitivity to surface
% sources and accounting for mixing histories", submitted, JPO.
% TMI = Total Matrix Intercomparison, Gebbie & Huybers, JPO, 2010.

%% Load the TMI tendency matrix, At, 
%  which satisfies dc/dt = At*c, where c is a global tracer distribution.
load L_4deg_2012

%% List the years for which the global tracer is output.
years = [0:1:100 110:10:1000 1025:25:2000 2100:100:5000];

yearsmid = (years(1:end-1)+years(2:end))./2;
NY = length(years)

%% Define the initial conditions.  
%  The surface layer is defined as all points in the mixed layer 
%  i.e, find(inmixlyr==1).
%  The time tendency in the surface layer defined to be zero.

  % Use the script 'make_initial_conditions',
% Or choose one of my pre-defined surface regions:
% 1) GLOBAL,2) ANT, 3 SUBANT, 4 NATL, 5 NPAC, 
% 6) TROP, 7. ARC, 8 MED, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
% 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
% 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP

% The next lines use the pre-defined regions. 
load c_all_4deg
wm = 1; % these are the pre-defined regions.
c0 = squeeze(c_all(:,wm));

% Pre-allocate arrays.
C = nan(NY,Nfield);
T = nan(NY,1);

% set the initial conditions.
C(1,:) = c0;

% options for the ODE solver.
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:Nfield,'Jacobian',At);

NN = (NY-1)./2;
for nn = 1:NN
  nn
  ind = nn*2-1:nn*2+1
  yrs = years(ind)
  Cin  = sq(C(ind(1),:));
  tic
  [Ttmp,Cout] = ode15s(@(t,x) At*x,yrs,Cin,options);
  toc
  C(ind,:) = sq(Cout);
  T(ind) = sq(Ttmp);
  %save ttdsummary C T % save intermediate values if necessary
end
save transient_output C T

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













