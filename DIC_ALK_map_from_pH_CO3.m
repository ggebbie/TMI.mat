%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        Find the distribution of a two tracers given:            %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC1 (tracer 1),   %
%           and dC2 (tracer2)
%           that best fits THE SPARSE observations, Cobs,         %
%           in tracer3 that is a nonlinear function of tracers 1,2 %
% Mathematically, minimize J = (E[C1,C2]-C3obs)^T W (E[C1,C2]-C3obs) subject to    %
%                         A C1 = dC1 + Gamma u1          AND       %
%                         A C2 = dC2 + Gamma u2          AND       %       
%  where u1 and u2 are the estimated changes in surface concentration.    % 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load a circulation field.
 load A_4deg_2012

LONv = LON(i); % grid longitudes in vector form.
LATv = LAT(j);

% get the longitude, latitude, depth, value of the
% observations, and their expected error.
% (lon_obs, lat_obs, depth_obs, y, err_obs)
lon_CO3_obs = [160 180 330 ];
lat_CO3_obs = [0   40   0  ];
depth_CO3_obs=[2e3 2e3 2e3 ];
%CO3_obs = [160 180 220]'; % picked out of thin air.
CO3_obs = [];
err_CO3_obs = 20;
if ~isempty(CO3_obs)
  WCO3  = speye(length(CO3_obs))./err_CO3_obs.^2;
  ECO3 = Emaker(lon_CO3_obs,lat_CO3_obs,depth_CO3_obs,LONv,LATv,k,DEPTH,6);
else
  WCO3 = [];
  ECO3 = [];
end

lon_pH_obs = [160 180 330 ];
lat_pH_obs = [0   40   0  ];
depth_pH_obs=[2e3 2e3 2e3 ];
pH_obs = [7.9 7.8 7.7]'; % picked out of thin air.
                         %pH_obs = [];
err_pH_obs = 0.1;
if ~isempty(pH_obs)
  WpH = speye(length(pH_obs))./err_pH_obs.^2;
  EpH = Emaker(lon_pH_obs,lat_pH_obs,depth_pH_obs,LONv,LATv,k,DEPTH,6);
else
  WpH = [];
  EpH = [];
end
W = blkdiag(WCO3,WpH);

y = [CO3_obs; pH_obs];

% set up the control variables.
Nsfc = sum(k==1);
isfc = find(k==1);
Nfield = length(k);
Gamma = sparse(isfc,1:Nsfc,ones(Nsfc,1),Nfield,Nsfc);

% what is the uncertainty in the surface boundary condition of DIC?
err_DIC_sfc = 200;
diagonal = false; % or set to false.
if diagonal 
  SDIC     = 1./err_DIC_sfc.^2; 
else
  % if diagonal==false
  % impose spatial smoothing in surface b.c.
  lengthscale = 4; % horizontal lengthscale in units of gridcells
  factor = 0.126.*(1/lengthscale)^2 ; 
  load Del2_4deg.mat
  SDIC = factor.* sparse(1:Nsfc,1:Nsfc,1./err_DIC_sfc.^2);
  SDIC = SDIC + Del2'*(lengthscale^4.*SDIC*Del2);
end
  
% what is the uncertainty in the surface boundary condition of TALK?
err_TALK_sfc = 200;
diagonal = false; % or set to false.
if diagonal 
  STALK     = 1./err_TALK_sfc.^2; 
else
  % if diagonal==false
  % impose spatial smoothing in surface b.c.
  lengthscale = 4; % horizontal lengthscale in units of gridcells
  factor = 0.126.*(1/lengthscale)^2 ; 
  load Del2_4deg.mat
  STALK = factor.* sparse(1:Nsfc,1:Nsfc,1./err_TALK_sfc.^2);
  STALK = STALK + Del2'*(lengthscale^4.*STALK*Del2);
end
S = blkdiag(SDIC,STALK);

% Get the first guess of surface boundary conditions and interior
% sources and sinks.
inmixlyr = find(diag(A)==1);
Nfield = length(i);
notmixlyr = (1:Nfield)'; notmixlyr(inmixlyr) = [];
 

% d0DIC= first guess of r.h.s. for DIC
dDIC0 = zeros(Nfield,1);   
dDIC0(k==1) = 2000;
dDIC0(notmixlyr) = -131.*dP; % from modern-day 
% cDIC0 = first guess DIC field
cDIC0 = A\dDIC0;

% d0TALK= first guess of r.h.s. for TALK
dTALK0 = zeros(Nfield,1);   
dTALK0(k==1) = 2200;
dTALK0(notmixlyr) = -30.*dP; % from modern-day 
% cTALK0 = first guess TALK field
cTALK0 = A\dTALK0;


%% Get y0,  first guess misfit.
% Get first guess CO3.
%cCO30 = CO3_from_DIC_ALK(cDIC0,cTALK0,Sal,TempC,P,Si,PO4,options)
pressure = DEPTH(k);
% Would be good to put actual T,S, etc. in next line.
addpath ~/mcode/toolbox/CO2SYS/

if ~isempty(CO3_obs)
  yCO3tilde0 = CO3_from_DIC_ALK(ECO3*cDIC0,ECO3*cTALK0,35,10,ECO3* ...
                                pressure,50,1,[]);
else
    yCO3tilde0 = [];
end
if ~isempty(pH_obs)
  ypHtilde0   = pH_from_DIC_ALK(EpH*cDIC0,EpH*cTALK0,35,10,EpH*pressure,50,1,[]);
else
  ypHtilde0 = [];
end

ytilde0 = [yCO3tilde0;ypHtilde0];
y0 = ytilde0-y; 

% check the cost function
Jobs0 = y0'*W*y0

cDIC = cDIC0;
cTALK = cTALK0;

J0 = Jobs0;
J = J0;

%% Nonlinear problem, so it may help to iterate a few times.
while J > 2

  if ~isempty(CO3_obs)  
    [EDIC,ETALK] = get_E_CO3(ECO3,cDIC,cTALK,pressure);

    %% Get F matrix from LaTex notes.
    FDIC_CO3  = get_F(A,Gamma,EDIC);
    FTALK_CO3 = get_F(A,Gamma,ETALK);
 
    FCO3 = [FDIC_CO3 FTALK_CO3];
  else
    FCO3 = [];
  end
  
  if ~isempty(pH_obs)  
    [EDIC,ETALK] = get_E_pH(EpH,cDIC,cTALK,pressure);

    %% Get F matrix from LaTex notes.
    FDIC_pH  = get_F(A,Gamma,EDIC);
    FTALK_pH = get_F(A,Gamma,ETALK);
 
    FpH = [FDIC_pH FTALK_pH];
  else
    FpH = [];
  end
  
  F = [FCO3; FpH];
  
  % get_hessian
  H = F'*(W*F) + S;

  f = F'*(W*y0);

  % direct solve.
  u = -H\f;

  %% calculate updated property fields.

  % DIC
  
  % alpha = 0.5; to eliminate a computational mode in gradient
  % descent.
  alpha = 1;
  dDIC = dDIC0 + alpha.*Gamma*u(1:Nsfc);
  cDIC = A\dDIC;

  % TALK
  dTALK = dTALK0 + alpha.*Gamma*u(Nsfc+1:end);
  cTALK = A\dTALK;

  if ~isempty(CO3_obs)
    yCO3tilde = CO3_from_DIC_ALK(ECO3*cDIC,ECO3*cTALK,35,10,ECO3*pressure,50,1,[]);
  else
      yCO3tilde = [];
  end
  if ~isempty(pH_obs)
    ypHtilde  = pH_from_DIC_ALK(EpH*cDIC,EpH*cTALK,35,10,EpH*pressure,50,1,[]);
  else
      ypHtilde = [];
  end
  
  ytilde = [yCO3tilde;ypHtilde];
  misfit = ytilde-y; 

  % check the cost function
  Jobs = misfit'*W*misfit

  Jcontrol = u'*S*u
  J = Jobs+Jcontrol
end


% how much did the data points reduce the error (globally).
% $$$ sqrt(sum((c-Tobs).^2)/Nfield)
% $$$ sqrt(sum((c0-Tobs).^2)/Nfield)
% $$$ 
% $$$ % how much did the data points reduce the error (at data points).
% $$$ Nobs = length(y);
% $$$ sqrt(sum((E*c-y).^2)/Nobs)
% $$$ sqrt(sum((E*c0-y).^2)/Nobs)
