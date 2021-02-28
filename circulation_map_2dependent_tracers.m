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
% load A_4deg_2012

% get the longitude, latitude, depth, value of the
% observations, and their expected error.
% (lon_obs, lat_obs, depth_obs, y, err_obs)
lon_obs = [160 180 330 ];
lat_obs = [0   40   0  ];
depth_obs=[2e3 2e3 2e3 ];


LONv = LON(i); % grid longitudes in vector form.
LATv = LAT(j);
E = Emaker(lon_obs,lat_obs,depth_obs,LONv,LATv,k,DEPTH,6);

% CO32- observations.
y = [160 180 220]'; % picked out of thin air.
err_CO3_obs = 20;

W  = 1./err_CO3_obs.^2;

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
cCO30 = CO3_from_DIC_ALK(E*cDIC0,E*cTALK0,35,10,E*pressure,50,1,[]);
%cCO30 = CO3_from_DIC_ALK(E*cDIC0,E*cTALK0,35.*ones(Nfield,1),10.*ones(Nfield,1),pressure,50.*ones(Nfield,1),ones(Nfield,1),[]);
y0 = cCO30-y; 

% check the cost function
Jobs0 = y0'*W*y0

% CO3 vs. TALK/DIC function.
aa = 37.0;
bb = 0.406;
cc = 0.00057

% $$$ % TALK minus DIC.
% $$$ AmD = cTALK0-cDIC0;
% $$$ 
% $$$ % Get EDIC. - linearizes around the first guess.
% $$$ e1 = speye(Nfield)*(-bb -(2.*cc).*AmD); 
% $$$ Etmp = sparse(1:Nfield,1:Nfield,e1);
% $$$ EDIC = E*Etmp;
% $$$ 
% $$$ % Get ETALK. - again linearized.
% $$$ e1 = speye(Nfield)*(bb +(2.*cc).*AmD); 
% $$$ Etmp = sparse(1:Nfield,1:Nfield,e1);
% $$$ ETALK = E*Etmp;

  % My functional form of CO3 vs (TALK-DIC) didn't work for the
  % sensitivities.
  % Instead try an impulse response method.
  [iii,jjj,sss]=find(E);

  sensDIC = zeros(length(ipert),1);
  ipert = unique(jjj);
  for ip = 1:length(ipert)
      delta = 1e-2;
    tmp0 =  CO3_from_DIC_ALK(cDIC0(ipert(ip)),cTALK0(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    tmp1 =  CO3_from_DIC_ALK(cDIC0(ipert(ip))+delta,cTALK0(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    sensDIC(ip) = (tmp1-tmp0)./delta;
    %senstest(ip) = -bb-2.*cc.*cCO3(ipert(ip));
  end
  Etmp = sparse(ipert,ipert,sensDIC,Nfield,Nfield);
  EDIC = E*Etmp;
   
  sensALK = zeros(length(ipert),1);
  for ip = 1:length(ipert)
      delta = 1e-2;
    tmp0 =  CO3_from_DIC_ALK(cDIC0(ipert(ip)),cTALK0(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    tmp1 =  CO3_from_DIC_ALK(cDIC0(ipert(ip)),cTALK0(ipert(ip))+delta,35,10,pressure(ipert(ip)),50,1,[]); 
    sensALK(ip) = (tmp1-tmp0)./delta;
  end
  Etmp = sparse(ipert,ipert,sensALK,Nfield,Nfield);
  ETALK = E*Etmp;


%% Get F matrix from LaTex notes.

% FDIC
if size(Gamma,2) <= size(EDIC,1)
  tic;
  tmp0 = A\Gamma;
  FDIC =  EDIC*tmp0;
  clear tmp0
  toc;
else
  tic;
  tmp0 = A'\EDIC';
  FDIC = transpose(Gamma'*tmp0);
  toc;
end

% FTALK
if size(Gamma,2) <= size(ETALK,1)
  tic;
  tmp0 = A\Gamma;
  FTALK =  ETALK*tmp0;
  clear tmp0
  toc;
else
  tic;
  tmp0 = A'\ETALK';
  FTALK = transpose(Gamma'*tmp0);
  toc;
end

%
F = [FDIC FTALK];
S = blkdiag(SDIC,STALK);

% get_hessian
H = F'*(W*F) + S;

f = F'*(W*y0);

% direct solve.
u =- H\f;

%% calculate updated property fields.

% DIC
dDIC = dDIC0 + Gamma*u(1:Nsfc);
cDIC = A\dDIC;

% TALK
dTALK = dTALK0 + Gamma*u(Nsfc+1:end);
cTALK = A\dTALK;

cCO3 = CO3_from_DIC_ALK(E*cDIC,E*cTALK,35,10,E*pressure,50,1,[]);


misfit = y-cCO3;
Jobs = misfit'*W*misfit;
Jcontrol = u'*S*u
J = Jobs+Jcontrol


%%  STRANGE: TRIED TO ITERATE BUT IT DIDN'T HELP. 
%%  wOULD RECOMMEND TO STOP HERE. 
%% IDEA: USE TARANTOLA AND VALETTE'S METHOD TO ITERATE, RATHER THAN
%% THE ONE i WROTE UP.
%% Nonlinear problem, so it may help to iterate a few times.
while J > 4
  % d0DIC= first guess of r.h.s. for DIC

  % update the first guess. 
  % mathematically not quite right because control variables have a
  % different reference point which is not taken into account in
  % the cost function.
  dDIC0 = dDIC;
  dTALK0 = dTALK;
  
  cDIC0 = cDIC:
  cTALK0 = cTALK;
 
  cCO30 = cCO3;

  % TALK minus DIC.
% $$$   AmD = cTALK0-cDIC0;
% $$$ 
% $$$   % Get EDIC. - linearizes around the first guess.
% $$$   e1 = speye(Nfield)*(-bb -(2.*cc).*AmD); 
% $$$   Etmp = sparse(1:Nfield,1:Nfield,e1);
% $$$   EDIC = E*Etmp;
% $$$ 
% $$$   % Get ETALK. - again linearized.
% $$$   e1 = speye(Nfield)*(bb +(2.*cc).*AmD); 
% $$$   Etmp = sparse(1:Nfield,1:Nfield,e1);
% $$$   ETALK = E*Etmp;

  % My functional form of CO3 vs (TALK-DIC) didn't work for the
  % sensitivities.
  % Instead try an impulse response method.
  [iii,jjj,sss]=find(E);
  
  ipert = unique(jjj);
  for ip = 1:length(ipert)
    delta = 1e-2;
    tmp0 =  CO3_from_DIC_ALK(cDIC(ipert(ip)),cTALK(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    tmp1 =  CO3_from_DIC_ALK(cDIC(ipert(ip))+delta,cTALK(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    sensDIC(ip) = (tmp1-tmp0)./delta;
    %senstest(ip) = -bb-2.*cc.*cCO3(ipert(ip));
  end
  Etmp = sparse(ipert,ipert,sensDIC,Nfield,Nfield);
  EDIC = E*Etmp;
   
  for ip = 1:length(ipert)
      delta = 1e-2;
    tmp0 =  CO3_from_DIC_ALK(cDIC(ipert(ip)),cTALK(ipert(ip)),35,10,pressure(ipert(ip)),50,1,[]); 
    tmp1 =  CO3_from_DIC_ALK(cDIC(ipert(ip)),cTALK(ipert(ip))+delta,35,10,pressure(ipert(ip)),50,1,[]); 
    sensALK(ip) = (tmp1-tmp0)./delta;
  end
  Etmp = sparse(ipert,ipert,sensALK,Nfield,Nfield);
  ETALK = E*Etmp;
    
  %% Get F matrix from LaTex notes.

  % FDIC
  if size(Gamma,2) <= size(EDIC,1)
    tic;
    tmp0 = A\Gamma;
    FDIC =  EDIC*tmp0;
    clear tmp0
    toc;
  else
    tic;
    tmp0 = A'\EDIC';
    FDIC = transpose(Gamma'*tmp0);
    toc;
  end

  % FTALK
  if size(Gamma,2) <= size(ETALK,1)
    tic;
    tmp0 = A\Gamma;
    FTALK =  ETALK*tmp0;
    clear tmp0
    toc;
  else
    tic;
    tmp0 = A'\ETALK';
    FTALK = transpose(Gamma'*tmp0);
    toc;
  end
 
  %
  F = [FDIC FTALK];

  % get_hessian
  H = F'*(W*F) + S;

  f = F'*(W*y0);

  % direct solve.
  u = H\f;

  %% calculate updated property fields.

  % DIC
  dDIC = dDIC0 + Gamma*u(1:Nsfc);
  cDIC = A\dDIC;

  % TALK
  dTALK = dTALK0 + Gamma*u(Nsfc+1:end);
  cTALK = A\dTALK;

  cCO3 = CO3_from_DIC_ALK(cDIC,cTALK,35.*ones(Nfield,1),10.*ones(Nfield,1),pressure,50.*ones(Nfield,1),ones(Nfield,1),[]);

  misfit = E*cCO3-y;
  Jobs = misfit'*W*misfit;
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
