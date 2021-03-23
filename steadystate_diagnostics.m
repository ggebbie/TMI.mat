%% steadystate_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
%  example 1: track surface waters into interior
%  example 2: Find the volume that has originated from surface
%  example 3: Find surface origin of an interior box
%  example 4: Reconstruct steady-state tracer field.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup the TMI data.

% Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

% Or, download manually at: https://docs.google.com/uc?export=download&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD and https://docs.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm

% Download TMI_4x4deg_data.tar.gz from Google Drive
! wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1L5i5eQ0QCrrqPKGoAxuB8X4_CltvQBMD' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm" -O TMI_4x4deg_data.tar.gz && rm -rf /tmp/cookies.txt

% Download TMI_2x2deg_data.tar.gz from Google Drive
! wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1xAkrTNybqoAKtFMuJ9XU9z9KZwnLenzm" -O TMI_2x2deg_data.tar.gz && rm -rf /tmp/cookies.txt


!tar xvzf TMI_2x2deg_data.tar.gz 
!rm TMI_2x2deg_data.tar.gz

!tar xvzf TMI_4x4deg_data.tar.gz 
!rm TMI_4x4deg_data.tar.gz

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
Nint = length(inotmixlyr)  % number of interior points

% if remineralized phosphate field is bigger than number of interior points, trim it for consistency.
if length(dP) > Nint
    dP = dP(inotmixlyr);
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
 if isempty(isec) 
     isec = find(LON==329);
 end
     
 figure(101)
 contourf(LAT,-DEPTH,squeeze(C(:,:,isec)),0:0.05:1) % a sample plot at 22 W.
 xlabel('latitude [deg N]')
 ylabel('depth [m]')
 
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
 if isempty(isec) 
     isec = find(LON==329);
 end
 figure(102)
 contourf(LAT,-DEPTH,squeeze(DPfld(:,:,isec)),0:0.05:1) % a sample plot at 22 W.
 xlabel('latitude [deg N]')
 ylabel('depth [m]')
 
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
 disty = 111.12 %km
 
 %for ny = 1:length(LAT)
 %  distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km');
 %end
 distx = transpose(dx.*disty.*cosd(LAT)); % formula instead of sw_dist
 
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
 xlabel('longitude')
 ylabel('latitude')
 
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
 N = size(A,1);
 
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
 contourf(LON,LAT,Vtot,[0 exp(-9:.5:-5)]) % quick plot

 xlabel('latitude')
 ylabel('longitude')
   
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
switch TMIversion
  case 'GH2010_2x2deg'
    load tracerobs_2deg_33lev_woce
  case 'GH2010_4x4deg'
    load tracerobs_4deg_33lev_woce
  case 'G2012_4x4deg'
    load tracerobs_4deg_33lev_woce
  case 'G2014_4x4deg'
    load tracerobs_4deg_33lev_woce
  otherwise
    disp('option not available')
end
  
  dT = Tobs;      % set rhs vector to the observations.
  dT(k>1) = 0;    % no internal sinks or sources.

  % a first guess: observed surface boundary conditions are perfect.  
  Tmod1 =  Q * (U \ (L \ (P * (R \ dT))));
  Tmod1_field = vector_to_field(Tmod1,i,j,k);
 
  % get cost function (J) for first guess temperature field.
  % here the data-model misfit is weighted by the expected error. 
  Nfield = length(Tobs);

  invWT{1} = sparse(1:Nfield,1:Nfield,1./Terr.^2,Nfield,Nfield); % diagonal matrix.
  J1 = (Tmod1-Tobs)'*invWT*(Tmod1-Tobs)./Nfield % optimized J1 = 1
 
  % Get better fit to observations by modifying the surface boundary 
  % condition within its uncertainty.
  % mathematical form: Find uT such that:
  % J = (Tmod-Tobs)'*inv(WT)*(Tmod-Tobs) is minimized 
  % subject to:  A*Tmod = dT + Gamma * uT.
  Nsfc = sum(k==1);
  Gamma = sparse(Nfield,Nsfc);
  nu = 0;
  for nv = 1:Nfield
    if k(nv)==1
      nu = nu+1;
      Gamma(nv,nu) = 1;
    end
  end
  u0 = zeros(Nsfc,1); % first guess of change to surface boundary
                      % conditions.
  lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
  ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

  %% 3 methods: 1) Constrained minimization with Lagrange multipliers 
  %(but without inequality constraints), 2) quadratic programming using full Hessian,
  % 3) quadratic programming using functional form of Hessian, 
  % 4) Unconstrained minimization according to Gebbie et al. 2015 (QSR). 

  % Here we proceed with method #1.
  options = optimset('Algorithm','interior-point','Display','iter', ...
                  'GradObj','on','LargeScale','on');
  % also can try 'Algorithm','trust-region-reflective'
  noncons = 0;
  isfc = find(kt==1);
  uTtilde= fmincon(@(x)objfun(x,A,invWT,Tobs,isfc,inotmixlyr,noncons),Tobs(isfc),[],[],[],[],lbT,ubT,[],options);
  J2 = objfun(uTtilde,A,invWT,Tobs,isfc,inotmixlyr,noncons) % expect J2 ~ 1, J2<J1
  d2 = zeros(Nfield,1); d2(isfc) = uTtilde;
  Tmod2 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
