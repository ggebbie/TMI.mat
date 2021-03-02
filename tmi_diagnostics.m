%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
%  example 1: track surface waters into interior
%  example 2: Find the volume that has originated from surface
%  example 3: Find surface origin of an interior box
%  example 4: Reconstruct steady-state tracer field.
%  example 5: Construct the pathways matrix from observations.
%  example 6: Follow the adjoint transport pathways.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup the TMI data.

% Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB.
%!wget https://drive.google.com/uc?export=download&id=1B03OHj82XT7A_ip48WtLmDF01zvQiz4v
%!wget https://drive.google.com/uc?export=download&id=1nFDSINbst84pK68aWwRGBVfYZkfN1mUR
%!wget https://drive.google.com/uc?export=download&id=1nFDSINbst84pK68aWwRGBVfYZkfN1mUR
%!wget https://drive.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm

!wget  --no-check-certificate 'https://drive.google.com/uc?export=download&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD' -O TMI_2x2deg_data.tar.gz

!wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm' -O TMI_4x4deg_data.tar.gz


!tar xvzf tmi_data.tar.gz 
!rm tmi_data.tar.gz

% make sure matlab paths are set
codedir = pwd
datadir = [codedir,'/data']

addpath(codedir) 
addpath(datadir)

%% Choose the TMI version of interest.
% For steady-state water-mass matrix, there are four choices.
TMIproducts = {'GH2010_2x2deg','GH2010_4x4deg','G2012_4x4deg','G2014_4x4deg'}

% choose the case
TMIno = 2; % select 1, 2, 3, 4
TMIversion = TMIproducts{TMIno}

switch TMIversion
  case 'GH2010_2x2deg'
    load A_2deg_2010
  case 'GH2010_4x4deg'
    load A_4deg_2010
  case 'G2012_4x4deg'
    load A_4deg_2012
  case 'G2014_4x4deg'
    load A_4deg_2014
  otherwise
    disp('option not available')
end

%% variables in A*.mat %%%
%
% A: pathways matrix
% i,j,k (or it,jt,kt): vectors that contain the spatial coordinates 
%        consistent with the A matrix.
% LAT,LON,DEPTH: the latitude, longitude and depth of the 3-d grid. 
% dP = amount of remineralized phosphate in each box (dimension either equal to all grid cells or only the grid cells below the mixed-layer depth)

%% make different versions of TMI consistent

% standardize the gridcell coordinates
if exist('it')
    i = it; j = jt; k = kt;
end

isfc = find(k==1); % surface gridcell indices
Nsfc = length(isfc); % number of surface points.

% Determine where the mixed layer is located in the vector of coordinates.
inmixlyr = find(diag(A)==1);
N = length(i) % total number of gridcells
Nmix = length(inmixlyr) % number of mixed layer points
inotmixlyr = (1:N)'; inotmixlyr(inmixlyr) = [];
Nint = length(notmixlyr)  % number of interior points

% if remineralized phosphate field is bigger than number of interior points, trim it for consistency.
if length(dP) > Nint
    dP = dP(iinotmixlyr);
end

%% Most examples are computationally more efficient with the output of a LU decomposition.
% Use an LU decomposition so that the inverse does not have
% to be stored explicitly.
% ALL EXAMPLES USE THIS CALCULATION.
[L, U, P, Q, R] = lu (A) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: Track the pathways of a user-defined water mass.         %
% Steps: (a) define the water mass 1). by a rectangular surface patch %
%            dyed with passive tracer concentration of 1,             % 
%            or 2. load a pre-defined surface patch in d_all.mat.     %
%        (b) propagate the dye with the matrix A, with the result     % 
%            being the fraction of water originating from the         %
%            surface region.                                          %
% See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the surface patch by the bounding latitude and longitude.
 lat_lo = 50; % 50 N, for example.
 lat_hi = 60;
 
 lon_lo = -50;
 lon_hi = 0;
 
 if lon_lo<=0 ;     lon_lo = lon_lo + 360; end
 if lon_hi<=0 ;     lon_hi = lon_hi + 360; end
 
 loc = find( LAT(j) < lat_hi & LAT(j)> lat_lo & ...
             LON(i) > lon_lo & LON(i)<lon_hi & k==1); 
 %potential wraparound problems in longitude are not accounted for
 %in the "loc" statement.
 
 N = length(i);
 d = zeros(N,1);
 d(loc) = 1; % a vector with ones in the N. Atl., zero elsewhere.
 
 % Also, can read in predefined surface patches in the file
 % d_all.mat, where the surface patches are defined as
 % oceanographically-relevant regions: 
 % 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
 % 6) ARC, 7 MED, 8) TROP, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
 % 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
 % 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP
 %
 % e.g.,
 % load d_all_4deg; d = d_all(:,3); %NATL.
 
 % to get quantity of dyed water throughout ocean:
 c = Q * (U \ (L \ (P * (R \ d))));  
 C = vector_to_field(c,i,j,k);  % g is the 3-d field of water-mass concentration.
 
 lon_section = 330;
 isec = find(LON==330);
 figure(101)
 contourf(LAT,-DEPTH,squeeze(C(:,:,isec)),0:0.05:1) % a sample plot at 22 W.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1B: Determine the total amount of remineralized phosphate   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % dP is total amount of remineralized phosphate in boxes below the
 % mixed-layer.
 d = zeros(N,1);
 d(inotmixlyr) = -dP; % a vector with remineralized phosphate sources.
 
 % DP is total remineralized phosphate including upstream contributions.
 % to get quantity of dyed water throughout ocean:
 DP = Q * (U \ (L \ (P * (R \ d))));  
 DPfld = vector_to_field(DP,i,j,k);  % g is the 3-d field of water-mass concentration.
 
 lon_section = 330;
 isec = find(LON==330);
 figure(102)
 contourf(LAT,-DEPTH,squeeze(DPfld(:,:,isec)),0:0.05:1) % a sample plot at 22 W.
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.
%
% See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 NZ = max(k);
 NY = max(j);
 NX = max(i);
 
 % get vector of volumes.
 dx = diff(LON(1:2));
 dy = diff(LAT(1:2));

 % km for 1 degree of latitude, could use sw_dist to get exact value.
 % disty = sw_dist([0 dy],[0 0],'km')
 disty = 111.12 

 
 %for ny = 1:length(LAT)
 %  distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km');
 %end
 distx = transpose(dx.*disty.*cosd(LAT)) % formula instead of sw_dist
 
 distx = reshape(distx,length(distx),1); % get dimensions correct.
 area = 1e6.*disty(ones(NY,NX)).*distx(:,ones(NX,1)); 
 zface= (DEPTH(1:end-1)+DEPTH(2:end))./2;
 dz = ([zface(1) ; diff(zface); 500]);
 vol = permute(area(:,:,ones(NZ,1)),[3 1 2]) .*dz(:,ones(NY,1),ones(NX,1));
 v   = field_to_vector(vol,i,j,k);
 
 % effectively take inverse of transpose A matrix.
 vtot = R' \ (P' * (L' \ (U' \ (Q' * v))));  % inverse of transpose
                                             % = transpose of inverse
 
 Vtot = vector_to_field(vtot,i,j,k);
 Vtot = sq(Vtot(1,:,:))./area; % scale the volume by surface area.
 figure(103)
 contourf(LON,LAT,log10(Vtot),1:0.25:7) % units: effective thickness
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3: Find the surface origin of water for some interior box %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of a
% given interior box,
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.      
% Very similar mathematically to example 2.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 NZ = max(k);
 NY = max(j);
 NX = max(i);
 Nfield = size(A,1);
 
 % choose an interior location X (Xlon[lon], Xlat [lat], Xdepth [m depth]).
% -7.38N, 115.26E
 Xlon = 125.26; % deg E.
 Xlat = -6.38; 
 Xdepth = 2500; % meters.
 
 % Find the coordinate on the grid by linear interpolation. 
 LONtmp = [LON(1)-4 ; LON; LON(end)+4];  %% watch out for the
                                         %wraparound of the grid.
 iX = interp1(LONtmp,0:length(LON)+1,Xlon,'linear')
 jX = interp1(LAT,1:length(LAT),Xlat,'linear')
 kX = interp1(DEPTH,1:length(DEPTH),Xdepth,'linear')
 kX(Xdepth>5500) = 33; % if deeper than the bottom level, put it on
                       % the bottom level.

 % watch out for wraparound again.
 i2 = i;
 i2( i-iX > NX/2) = i2( i-iX > NX/2) - NX;
 i2( iX-i > NX/2) = i2( iX-i > NX/2) + NX;
  
 % find the closest gridpoints to the location of interest.
 [dist,loc] = sort((i2-iX).^2 + (j-jX).^2 + (k-kX).^2);
 dist = dist+0.1; % to eliminate singularity if point matches up
                  % exactly with the grid. 
 
 mask = zeros(N,1); 
 mask(loc(1:6)) =(1./dist(1:6))./sum(1./dist(1:6)); % mask adds up
                                                    % to 1.
 
 vtot = R' \ (P' * (L' \ (U' \ (Q' * mask)))); 
 Vtot = vector_to_field(vtot,i,j,k);
 
 Vtot = squeeze(Vtot(1,:,:)); % just keep the surface as that is
                              % the ultimate source region.
 figure(104)
 contourf(LON,LAT,Vtot,[0 .001 .005 .01 .05 .1 .5]) % one way to
                                                    % plot it.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4: Find the distribution of a tracer given:              %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC,               % 
%           that best fits observations, Cobs,                    %
%   and (c) inequality constraints on the tracer concentration.   %
%                                                                  %
% Mathematically, minimize J = (C-Cobs)^T W (C-Cobs) subject to    %
%                         AC = d + Gamma u                         %
%  where u is the estimated change in surface concentration.    % 
%
% See Supplementary Section 2, Gebbie & Huybers 2011.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  % use temperature as an example, with Tobs the observed temperature vector.
  % load temperature in vector form, along with associated error.
  % choose version of tracer observations
dataname = {'tracerobs_2deg_33lev_woce.mat','tracerobs_4deg_33lev_woce.mat','tracerobs_4deg_33lev_TMI.mat'}
datano = 2; % choose 1,2, or 3
if exist(dataname{datano})
    load(dataname{datano})
else
    disp('not an available choice')
end

% check that the choice is consistent with pathways matrix.
if length(Tobs) ~= N
    disp('inconsistent obs and pathways matrix')
end

dT = Tobs;      % set rhs vector to the observations.
dT(k>1) = 0;    % no internal sinks or sources.

% a first guess: observed surface boundary conditions are perfect.  
Tmod1 =  Q * (U \ (L \ (P * (R \ dT))));
Tmod1_field = vector_to_field(Tmod1,i,j,k);
 
% get cost function (J) for first guess temperature field.
% here the data-model misfit is weighted by the expected error. 
% N = length(Tobs);
invWT{1} = sparse(1:N,1:N,1./Terr.^2,N,N); % diagonal matrix.
J1 = (Tmod1-Tobs)'*invWT{1}*(Tmod1-Tobs)./N % optimized J1 = 1
 
% Get better fit to observations by modifying the surface boundary 
% condition within its uncertainty.
% mathematical form: Find uT such that:
% J = (Tmod-Tobs)'*inv(WT)*(Tmod-Tobs) is minimized 
% subject to:  A*Tmod = dT + Gamma * uT.
% Nsfc = sum(k==1);
Gamma = sparse(N,Nsfc);
nu = 0;
for nv = 1:N
    if k(nv)==1
        nu = nu+1;
        Gamma(nv,nu) = 1;
    end
end
u0 = zeros(Nsfc,1); % first guess of change to surface boundary
                      % conditions.
lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

%% choose a method below 

method = 1;
if method==1
    % 1) quadratic programming using full Hessian
    H1 = Q * (U \ (L \ (P * (R \ dT)))) -Tobs; 
    H2 = Q * (U \ (L \ (P * (R \ Gamma))));   % lengthy calculation step: get full
                                              % Hessian. Can be avoided by specifying
                                              % functional form of Hessian-vector product.

    tmp1 = invWT{1}*H2;
    HT = H2'*tmp1;
    fT = 2.*H1'*tmp1;
    uT = quadprog(HT,fT,[],[],[],[],lbT,[],u0); % solve quadratic
                                                % programming problem for
                                                % shifts (uT) to
                                                % sfc. temp.
    J1 = objfun(uTtest,A,invWT,Tobs,isfc,inotmixlyr,noncons) %
                                       
    d1 = zeros(N,1); d1(isfc) = uT;
    Tmod1 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
    
elseif method==3 % slower because doesn't use Hessian.
      %%  3) Constrained minimization (without inequality constraints),
      
      maxiters = 10; % increase for a better fit.
      options = optimset('Algorithm','interior-point','Display','iter', ...
                         'GradObj','on','LargeScale','on','maxiter',maxiters);
      %options = optimset('Algorithm','trust-region-reflective','Display','iter', ...
      %        'GradObj','on','LargeScale','on');
      %isfc = find(kt==1);

      noncons = 0;
      uTtest= fmincon(@(x)objfun(x,A,invWT,Tobs,isfc,inotmixlyr,noncons),Tobs(isfc),[],[],[],[],lbT,ubT,[],options);
      J2 = objfun(uTtest,A,invWT,Tobs,isfc,inotmixlyr,noncons) %
                                                         % expect J2 ~ 1, J2<J1
      d2 = zeros(N,1); d2(isfc) = uTtest;
      Tmod2 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
end
% Note: other methods are: 3) quadratic programming using functional form of Hessian, 
% 4) Unconstrained minimization. 
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4B: Same as example 4 but the observations are sparse and
% the inequality constraint is dropped.
%
%        Find the distribution of a tracer given:                 %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC,               % 
%           that best fits THE SPARSE observations, Cobs,         %
%                                                                  %
% Mathematically, minimize J = (C-Cobs)^T W (C-Cobs) subject to    %
%                         AC = d + Gamma u                         %
%  where u is the estimated change in surface concentration.    % 
%
% See Supplementary Section 2, Gebbie & Huybers 2011.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load a circulation field.
% load A_4deg_2012

% get the longitude, latitude, depth, value of the
% observations, and their expected error.
% (lon_obs, lat_obs, depth_obs, y, err_obs)
lon_obs = [160 180 200 220];
lat_obs = [0   0   0   0 ];
depth_obs=[2e3 2e3 2e3 2e3];


LONv = LON(i); % grid longitudes in vector form.
LATv = LAT(j);
E = Emaker(lon_obs,lat_obs,depth_obs,LONv,LATv,k,DEPTH,6);

% For this example, sample the observations perfectly from existing
% temperature observations.
load tracerobs_4deg_33lev_woce Tobs Terr
y = E*Tobs;
err_obs = 0.1;

W  = 1./err_obs.^2;

% set up the control variables.
Nsfc = sum(k==1);
isfc = find(k==1);
N = length(k);
Gamma = sparse(isfc,1:Nsfc,ones(Nsfc,1),N,Nsfc);

% what is the uncertainty in the surface boundary condition?
err_d = 1;
diagonal = false; % or set to false.
if diagonal 
  S     = 1./err_d.^2; 
else
  % if diagonal==false
  % impose spatial smoothing in surface b.c.
  lengthscale = 4; % horizontal lengthscale in units of gridcells
  factor = 0.126.*(1/lengthscale)^2 ; 
  load Del2_4deg.mat
  S = factor.* sparse(1:Nsfc,1:Nsfc,1./err_d.^2);
  S = S + Del2'*(lengthscale^4.*S*Del2);
end
  
% Get the first guess of surface boundary conditions and interior
% sources and sinks.
% d= first guess of r.h.s.
d0 = Tobs;   % truth is imposed as first-guess boundary condition.
d0(k>1) = 0; % conservative in interior

% c0 = first guess property field
c0 = A\d0;

% r = first guess misfit.
r = E*c0-y;

% get transpose(f)
f1 =   W*(E'*r);
f2 =  A'\f1;
f  = Gamma'*f2;
clear fT1 fT2

% 2 methods to solve: iterative or direct
iterative = false;
if iterative % method 1
  % put in the form of a least-squares problem. Hu = -f;
  u = lsqr(@(x,tflag)afun(x,Gamma,A,E,W,S),-f,1e-6,20);
else % method 2
  tic; H = get_hessian(Gamma,A,E,W,S); toc;
  u = -H\f;
end

% update the boundary conditions with the least squares solution.
d = d0 + Gamma*u;

% solve for the three dimensional property field.
c = A\d;
cfld = vector_to_field(c,i,j,k);

% how much did the data points reduce the error (globally).
sqrt(sum((c-Tobs).^2)/N)
sqrt(sum((c0-Tobs).^2)/N)

% how much did the data points reduce the error (at data points).
Nobs = length(y);
sqrt(sum((E*c-y).^2)/Nobs)
sqrt(sum((E*c0-y).^2)/Nobs)
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4C: Same as example 4B but the tracer is nonconservative.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load a circulation field.
% load A_4deg_2012

% get the longitude, latitude, depth, value of the
% observations, and their expected error.
% (lon_obs, lat_obs, depth_obs, y, err_obs)
lon_obs = [160 180 200 220];
lat_obs = [0   0   0   0 ];
depth_obs=[2e3 2e3 2e3 2e3];

LONv = LON(i); % grid longitudes in vector form.
LATv = LAT(j);
E = Emaker(lon_obs,lat_obs,depth_obs,LONv,LATv,k,DEPTH,6);

% For this example, sample the observations perfectly from existing
% temperature observations.
load tracerobs_4deg_33lev_woce Pobs 
y = E*Pobs;
err_obs = 0.1;

W  = 1./err_obs.^2;

% set up the control variables.

% what is the uncertainty in the first guess global field?
% for sake of argument, use the global observations as the first
% guess.
c0 = 1.0.*ones(N,1);
err_c0 = 3;
S     = 1./err_d.^2; % this step could be done better -- impose
                     % spatial smoothing in surface b.c.

% Get the first guess of surface boundary conditions and interior
% sources and sinks.
% d= first guess of r.h.s.

% diagnose the nonconservative effects in the dataset.
% surface boundary conditions.
q = A*Pobs;

% to be consistent: use same surface boundary condition from above.

% what is error in surface boundary conditions?
err_sfc = 5;

% what is the assumed error in the nonconservative effects?
err_noncons = 1e-3;

% turn into a weighting matrix.
qerr = zeros(N,1);
qerr(k==1) = err_sfc;
qerr(k>1) = err_noncons;
Q = sparse(1:N,1:N,1./qerr.^2);

ETWE = E'*(W*E);
ATQA = A'*(Q*A);

H = ETWE+ATQA+S*speye(N);

f = E'*(W*y) + A'*(Q*q) + S.*(c0);

c = H\f;

% update the boundary conditions with the least squares solution.
cfld = vector_to_field(c,i,j,k);

% how much did the data points reduce the error (globally).
sqrt(sum((c-Pobs).^2)/N)
sqrt(sum((c0-Pobs).^2)/N)

% how much did the data points reduce the error (at data points).
Nobs = length(y);
sqrt(sum((E*c-y).^2)/Nobs)
sqrt(sum((E*c0-y).^2)/Nobs)

% unlike Example 4, does not enforce non-negativity.

   
