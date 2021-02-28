%% ARCTIC CAP -- WORK WITH DAVID THORNALLEY. JULY 2012.
% EXP 0. CONTROL EXPERIMENT WITH C14. 
%        A. GET C14 DATA.
%        b. gET PATHWAYS AND TURNOVER TIME FOR EACH BOX.
%        c. sOLVE FOR MODERN-DAY C14 THAT BEST FITS DATA.
% EXP 1. PUT SOME "BOUNDARIES" AT THE BOTTOM OF THE ARCTIC OCEAN.
%       SET THEIR VALUE TO BE 400 PER MIL LESS THAN THEY ACTUALLY ARE.
%        A. KEEP OTHER BOUNDARIES AT MODERN-VALUES.
%        b. COMPUTE C14 EVERYWHERE.
%        c. LOOK AT DIFFERENCE COMPARED TO REFERENCE VALUE.
% EXP 2. CAP THE ARCTIC OCEAN.
%        a. CHANGE THE PATHWAYS TO REMOVE THE ARCTIC SURFACE AS 
%           a bOUNDARY POINT. MAKE WATER DIFFUSE AWAY EVENLY FROM
%           THE SURFACE. 
%         b. RE-RUN WITH THE MODERN-DAY BOUNDARY CONDITIONS. 
%        C. COMPARE TO MODERN-DAY REFERENCE CONDITIONS. SHOULD GIVE
%        SOME DIFFERENCE DUE TO EXTRA ISOLATION. 



%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
% load the master file: choose either woce2deg or woce4deg 
% (2 deg == 2 x 2 horiz resolution, 33 vertical levels)
% (4 deg == 4 x 4 horiz resolution, 33 vertical levels)
% (woce == trained on the WOCE Global Hydrographic Climatology,
%          Gouretski and Koltermann 2004)

load A_woce4deg  %% choose between the 4 deg and 2 deg MAT-files
 
% variables in A.mat %%%
%
% A: pathways matrix
% i,j,k: vectors that contain the spatial coordinates 
%        consistent with the A matrix.
% LAT,LON,DEPTH: the latitude, longitude and depth of the 3-d grid. 
% dP = amount of remineralized phosphate in each box.

%% Most examples are computationally more efficient with the output of a LU decomposition.
% Use an LU decomposition so that the inverse does not have
% to be stored explicitly.

load L_4deg_10feb2011

Nfield = length(it);
% Turn At into a radioactive decay matrix.
isfc = find(kt==1);
Tau = -1./diag(At);
Tau(abs(Tau)>1e5) = 0;
lambda = log(2)./5730;

T = sparse(1:Nfield,1:Nfield,Tau);

A14 = A - lambda*T;
[Lx, Ux, Px, Qx, Rx] = lu (A14) ;

load ~/mcode/c14/c14work_26july2011_cleaned.mat copt4b

dC = copt4b;
dC(kt>1) = 0;
 
corig = Qx * (Ux \ (Lx \ (Px * (Rx \ dC))));  
C = vector_to_field((corig-1).*1000,i,j,k);  % g is the 3-d field of water-mass concentration.

contourf(LAT,-DEPTH,squeeze(C(:,:,85))) % a sample plot at 22 W.
 
%% exp 1. arctic isolation.
load ~/stats_4deg_33lev_woce_4jan.mat atbottom

%% add new boundary points.
ibot = find(LAT(jt)>=67 & atbottom);
ibdy = [isfc; ibot];
 
dCiso = zeros(Nfield,1);
dCiso(ibdy) = c(ibdy);

A14iso = A14;
for ii= 1:length(ibot)
    ii
  A14iso(ibot(ii),:) = 0;
  A14iso(ibot(ii),ibot(ii)) = 1;
end
[Lx, Ux, Px, Qx, Rx] = lu (A14iso) ;

ciso = Qx * (Ux \ (Lx \ (Px * (Rx \ dCiso))));  
Ciso = vector_to_field((ciso-1).*1000,i,j,k);  % g is the 3-d field of water-mass conc 

std(ciso-corig)

%  Now lower C14 by 400 permil = 0.4
dCiso(ibot) = dCiso(ibot)-0.4;
ciso = Qx * (Ux \ (Lx \ (Px * (Rx \ dCiso))));  
Ciso = vector_to_field((ciso-1).*1000,i,j,k);  % g is the 3-d field
                                               % of water-mas

contourf(LAT,-DEPTH,squeeze(Ciso(:,:,85))) % a sample plot at 22 W.

Ciso = vector_to_field((ciso-1).*1000,i,j,k);  % g is the 3-d field
                                               % of water-mas

delCiso = Ciso-C;
contourf(LAT,-DEPTH,squeeze(delCiso(:,:,85))) % a sample plot at 22 W.

cbot = get_bottom_tracer(ciso-corig,atbottom,it,jt);
figure
[c,h]=contourf(LON,LAT,squeeze(1000.*cbot))

%% exp 2 -- cap the arctic.

load d_all_4deg
iarc = find(d_all(:,6) | d_all(:,12));

NX = 90;
A14cap = A14;
for ii = 1:length(iarc)
    nv = iarc(ii)
    % find nearest surface points and divide them up evenly.

    if it(nv)==90
      neigh = find((it==1 & jt ==jt(nv) & kt==1) | ...
                 (it==it(nv)-1 & jt ==jt(nv) & kt==1) |...
                 (it==it(nv) & jt==jt(nv)-1  & kt==1) |...
                 (it==it(nv) & jt==jt(nv)+1  & kt==1));
    elseif it(nv)==1
      neigh = find((it==it(nv)+1 & jt ==jt(nv) & kt==1) |...
                 (it==90 & jt ==jt(nv) & kt==1) |...
                 (it==it(nv) & jt==jt(nv)-1  & kt==1) |...
                 (it==it(nv) & jt==jt(nv)+1  & kt==1));
                 
    else
      neigh = find((it==it(nv)+1 & jt ==jt(nv) & kt==1) |...
                 (it==it(nv)-1 & jt ==jt(nv) & kt==1) |...
                 (it==it(nv) & jt==jt(nv)-1  & kt==1) |...
                 (it==it(nv) & jt==jt(nv)+1  & kt==1));
        
    end
        
     % find the closest gridpoints to the location of interest.
    A14cap(nv,:) = 0;     
    A14cap(nv,neigh) = 1./length(neigh);
    A14cap(nv,nv) = -1;
    
end
[Lx, Ux, Px, Qx, Rx] = lu (A14cap) ;
dCcap = copt4b;
dCcap(kt>1) = 0;
dCcap(iarc) = 0;
 
ccap = Qx * (Ux \ (Lx \ (Px * (Rx \ dCcap))));  
Ccap = vector_to_field((ccap-1).*1000,i,j,k);  % g is the 3-d field of water-mass concentration.

contourf(LAT,-DEPTH,squeeze(Ccap(:,:,85))) % a sample plot at 22 W.

%
delccap = ccap - corig;
delCcap = vector_to_field((delccap).*1000,i,j,k);  % g is the 3-d field of water-mass concentration.


%%% example 2: slow down the Arctic circulation by 4 (X) times.

Tau2 = Tau;
Tau2(LAT(jt)<67) = 0;
Tau2(LAT(jt)>=67) = 4.*Tau2(LAT(jt));
T2 = sparse(1:Nfield,1:Nfield,Tau2);

A14capslow = A14cap - lambda*T2;
[Lx, Ux, Px, Qx, Rx] = lu (A14capslow) ;
ccapslow = Qx * (Ux \ (Lx \ (Px * (Rx \ dCcap))));  
Ccapslow = vector_to_field((ccap-1).*1000,i,j,k);  % g is the 3-d field of water-mass concentration.

delccapslow = ccapslow - corig;
delCcapslow = vector_to_field((delccapslow).*1000,i,j,k);  % g is the 3-d field of water-mass conc
contourf(LAT,-DEPTH,squeeze(delCcapslow(:,:,85))) % a sample plot at 22 W.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: Track the pathways of a user-defined water mass.         %
% Steps: (a) define the water mass 1). by a rectangular surface patch %
%            dyed with passive tracer concentration of 1,             % 
%            or 2. load a pre-defined surface patch in d_all.mat.     %
%        (b) propagate the dye with the matrix A, with the result     % 
%            being the fraction of water originating from the         %
%            surface region.                                          %
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

 % to get quantity of dyed water throughout ocean:
 c = Q * (U \ (L \ (P * (R \ d))));  
 C = vector_to_field(c,i,j,k);  % g is the 3-d field of water-mass concentration.
 
 contourf(LAT,-DEPTH,squeeze(C(:,:,85))) % a sample plot at 22 W.
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 NZ = max(k);
 NY = max(j);
 NX = max(i);
 
 % get vector of volumes.
 dx = diff(LON(1:2));
 dy = diff(LAT(1:2));
 for ny = 1:length(LAT)
   distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km');
 end
 disty = sw_dist([0 dy],[0 0],'km')
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
 Vtot = squeeze(Vtot(1,:,:))./area; % scale the volume by surface area.
 contourf(LON,LAT,log(Vtot))

 
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
% -7.38, 115.26E
 Xlon = 125.26; % deg E.
 Xlat = -6.38; 
 Xdepth = 500; % meters.
 
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
 
 mask = zeros(Nfield,1); 
 mask(loc(1:6)) =(1./dist(1:6))./sum(1./dist(1:6)); % mask adds up
                                                    % to 1.
 
 vtot = R' \ (P' * (L' \ (U' \ (Q' * mask)))); 
 Vtot = vector_to_field(vtot,i,j,k);
 
 Vtot = squeeze(Vtot(1,:,:)); % just keep the surface as that is
                              % the ultimate source region.

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
%  where u is the estimated change in surface concentration.       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  % use temperature as an example, with Tobs the observed temperature vector.
  % load temperature in vector form, along with associated error.
  load tracerobs_4deg_33lev_woce  % could also choose 2deg

  
  dT = Tobs;      % set rhs vector to the observations.
  dT(k>1) = 0;    % no internal sinks or sources.

  % a first guess: observed surface boundary conditions are perfect.  
  Tmod1 =  Q * (U \ (L \ (P * (R \ dT))));
  Tmod1_field = vector_to_field(Tmod1,i,j,k);
 
  % get cost function (J) for first guess temperature field.
  % here the data-model misfit is weighted by the expected error. 
  Nfield = length(Tobs);
  Terr(ibotside) = Terr(ibotside).*50;
  %Terr(isfc) = 0;
  invWT = sparse(1:Nfield,1:Nfield,1./Terr.^2,Nfield,Nfield); % diagonal matrix.
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

  %% 2 methods: a) compute full Hessian, b) use functional form of
  %Hessian. (b) is preferred.
  method = 3;
  if method==1
    H1 = Q * (U \ (L \ (P * (R \ dT)))) -Tobs; 
    H2 = Q * (U \ (L \ (P * (R \ Gamma))));   % lengthy calculation step: get full
                                 % Hessian. Can be avoided by specifying
				 % functional form of Hessian-vector product.

    HT = H2'*invWT*H2;
    fT = 2.*H1'*invWT*H2;
    uT = quadprog(HT,fT,[],[],[],[],lbT,[],u0); % solve quadratic
                                                 % programming problem for
                                                 % shifts (uT) to
                                                 % sfc. temp.
    
  elseif method==2
    [fval, exitflag, output, uT] = hessianprob(A,Gamma,invWT,Tobs,u0,dT,lbT,ubT);
    Tmod2 =  Q * (U \ (L \ (P * (R \ (dT+Gamma*uT))))) ; % best estimate
    J2 = (Tmod2-Tobs)'*invWT*(Tmod2-Tobs)./Nfield % expect J2 ~ 1, J2<J1

  elseif method==3 % slower because doesn't use Hessian.
      %%
     options = optimset('Algorithm','interior-point','Display','iter', ...
                  'GradObj','on','LargeScale','on');
     %options = optimset('Algorithm','trust-region-reflective','Display','iter', ...
     %        'GradObj','on','LargeScale','on');
     noncons = 0;
     uTtest= fmincon(@(x)objfun(x,A,invWT,Tobs,isfc,iint,noncons),Tobs(isfc),[],[],[],[],lbT,ubT,[],options);
     J2 = objfun(uTtest,A,invWT,Tobs,isfc,iint,noncons) %
                      % expect J2 ~ 1, J2<J1
     d2 = zeros(Nfield,1); d2(isfc) = uTtest;
     Tmod2 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
  
  end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 5: How to update the TMI pathways and tracer distribution    %
%            to improve the fit to observations.                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% modify directions in Schlitzer 2007.

load ../woce_4deg/A0_woce_4deg_4jan
load ../woce_4deg/stats_4deg_33lev_woce_4jan.mat W*sfc *err
invWTsfc = WTsfc; clear WTsfc;
invWSsfc = WSsfc; clear WSsfc;
invWO18sfc = WO18sfc; clear WO18sfc;
invWPsfc = WPsfc; clear WPsfc;
invWNsfc = WNsfc; clear WNsfc;
invWOsfc = WOsfc; clear WOsfc;

%% Get all 6 tracers.
%% add something for ibotside?
Terr_nosfc = Terr; Terr_nosfc(isfc) = Inf;
invWTsm = sparse(1:Nfield,1:Nfield,1./Terr_nosfc.^2,Nfield,Nfield); % diagonal matrix.
invWTsm = invWTsm + Gamma*invWTsfc*Gamma'; % weighting matrix with
                                           % surface smoothing

Serr_nosfc = Serr; Serr_nosfc(isfc) = Inf;
invWSsm = sparse(1:Nfield,1:Nfield,1./Serr_nosfc.^2,Nfield,Nfield);
invWSsm = invWSsm + Gamma*invWSsfc*Gamma'; 

O18err_nosfc = O18err; O18err_nosfc(isfc) = Inf;
invWO18sm = sparse(1:Nfield,1:Nfield,1./O18err_nosfc.^2,Nfield,Nfield);
invWO18sm = invWO18sm + Gamma*invWO18sfc*Gamma';

Perr_nosfc = Perr; Perr_nosfc(isfc) = Inf;
invWPsm = sparse(1:Nfield,1:Nfield,1./Perr_nosfc.^2,Nfield,Nfield);
invWPsm = invWPsm + Gamma*invWPsfc*Gamma';

Nerr_nosfc = Nerr; Nerr_nosfc(isfc) = Inf;
invWNsm = sparse(1:Nfield,1:Nfield,1./Nerr_nosfc.^2,Nfield,Nfield);
invWNsm = invWNsm + Gamma*invWNsfc*Gamma';

Oerr_nosfc = Oerr; Oerr_nosfc(isfc) = Inf;
invWOsm = sparse(1:Nfield,1:Nfield,1./Oerr_nosfc.^2,Nfield,Nfield);
invWOsm = invWOsm + Gamma*invWOsfc*Gamma';

Wall{1} = invWTsm;
Wall{2} = invWSsm;
Wall{3} = invWO18sm;
Wall{4} = invWPsm;
Wall{5} = invWNsm;
Wall{6} = invWOsm;

Wpno{1} = invWPsm;
Wpno{2} = invWNsm;
Wpno{3} = invWOsm;

%% use objfun without changing paths -- but with 6 tracers and
%% nonconservative effects.



%%
m = A(iintdst+(iintsrc-1).*Nfield)';
x0 = [Tmod2(k==1); m];
Atemplate = A;
Atemplate(~inmixlyr,:) = 0;

D = sparse(iintdst,1:length(iintsrc),ones(Nm,1),Nfield,Nm);
D(inmixlyr,:) = [];

%%
%[Jtest1,gtest1] = objfun(x0,Atemplate,Wall,Tobs,isfc,iint,0, ...
%                         iintdst,iintsrc); 


% test case: take A as fixed. Compare J with optimized sfc + dp
% vs. without.

x0 = [Tobs(isfc);Sobs(isfc);O18obs(isfc);Pobs(isfc);Nobs(isfc); ...
      Oobs(isfc);dP(iint)];
cobs = [Tobs Sobs O18obs Pobs Nobs Oobs];
nclist = [0 0 0 1 15.5 -170];
[J0,g0] = objfun(x0,A,Wall,cobs,isfc,iint,nclist); 

delta = 1e-3; deltaloc = 50000;
x_path = x0;
x_path2 = x_path; x_path2(deltaloc) = x_path2(deltaloc) + delta;
[J2A,g2A] = objfun(x_path,A,Wall,cobs,isfc,iint,nclist);
[J2B,g2B] = objfun(x_path2,A,Wall,cobs,isfc,iint,nclist);
(g2B(deltaloc)+g2A(deltaloc))./2
(J2B-J2A)./delta

%% P,N, O optimization.
lb = [zeros(3*Nsfc,1); zeros(Nint,1)];
ub = [];

%ub = [35 + 0.*Tmod2(k==1); 1 + 0.*m];
%ub = [35 + 0.*Tmod2(k==1)];

while J2 > 1
   options = optimset('Algorithm','interior-point','Display','iter', ...
                   'GradObj','on','LargeScale','on','MaxIter',20, ...
                   'MaxFunEvals',40,'Hessian',{'lbfgs',4});
   %options = optimset('Algorithm','trust-region-reflective','Display','iter', ...
   %                   'GradObj','on','LargeScale','on','MaxIter',5, ...
   %                   'MaxFunEvals',10); % too much memory.

 [x_path,fval,exitflag,output]= fmincon(@(x)objfun(x,A, ...
     Wpno,sq(cobs(:,4:6)),isfc,iint,nclist(4:6)),x0(3*Nsfc+1:end), ...
     [],[],[],[],lb(3*Nsfc+1:end),ub,[],options);




   noncons = 0;
   [x_path,fval,exitflag,output]= fmincon(@(x)objfun(x,Atemplate,invWTsm,Tobs, ...
     isfc,iint,noncons,iintdst,iintsrc),x_path,[],[],[],[],lb,ub,[],options);

   %x_path= fmincon(@(x)objfun(x,Atemplate,invWT,Tobs, ...
   %  isfc,iint,noncons,iintdst,iintsrc),x0,[],[],[],[],lb,ub,[],options);

   [J2,g2] = objfun(x_path,Atemplate,invWTsm,Tobs,isfc,iint,noncons,iintdst,iintsrc); %
                      % expect J2 ~ 1, J2<J1

   % renormalize the "m" part of the state vector.
   diagsum = D*x_path(end-Nm+1:end);

   M = sparse(iintdst,iintsrc,x_path(end-Nm+1:end),Nfield,Nfield);
   Diagsum = sparse(1:Nint,1:Nint,1./diagsum,Nint,Nint);
   M(iint,:) = Diagsum*M(iint,:);
   
   m = M(iintdst+(iintsrc-1).*Nfield);
   x_path = [x_path(1:2806); m'];
   save work_1dec2011
end


%% check the gradient.                      
delta = 1e-3
deltaloc = 2000;

x_path2 = x_path; x_path2(deltaloc) = x_path2(deltaloc) + delta;
[J2A,g2A] = objfun(x_path,Atemplate,invWTsm,Tobs,isfc,iint,noncons,iintdst,iintsrc); 
[J2B,g2B] = objfun(x_path2,Atemplate,invWTsm,Tobs,isfc,iint, ...
                   noncons,iintdst,iintsrc); 
(g2B(deltaloc)+g2A(deltaloc))./2
(J2B-J2A)./delta

%%
load ../woce_4deg/grid_4deg_33lev_woce.mat

% Go through all interior locations to clean up indices.
iint = find(~inmixlyr);
imix = find(inmixlyr);

% clean up my destination and source info to not include mixed-layer
% destinations.
iintsrc(isnan(iintsrc))= [];
iintdst(isnan(iintdst))= [];
iintdst = transpose(iintdst);

for nm = 1:length(imix);
   
   nv = imix(nm)
   iintsrc(iintdst==nv)=[];
   iintdst(iintdst==nv)=[];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6: How to modify or reconstruct the TMI pathways,            %
%            i.e., the A matrix.                                       %
%                                                                      %
% Each row of A is determined by the solution of a non-negative        %
% least squares problem via the method of Lawson and Hanson.           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Pick an example box number "nv". Eventually must go through every box.  
% nv = 30000;

% If in mixed layer, then stop. The local inversion is not applicable.
if ~inmixlyr(nv)
  
  %% Find the neighboring gridcells.
  %  Account for longitudinal wraparound.
  it2 = it;
  it2( it-it(nv) > NX/2) = it2( it-it(nv) > NX/2) - NX;
  it2( it(nv)-it > NX/2) = it2( it(nv)-it > NX/2) + NX;
  boxsize =  1;

  %% six sides to a box.
  isource1=find( abs(it2-it2(nv))==boxsize & ...
              jt + boxsize > jt(nv)   & ...
              jt - boxsize < jt(nv)   & ...
              kt + boxsize > kt(nv)   & ...
              kt - boxsize < kt(nv)  );
  isource2=find( abs(jt-jt(nv))==boxsize & ...
              it2 + boxsize > it2(nv)   & ...
              it2 - boxsize < it2(nv)   & ...
              kt + boxsize > kt(nv)   & ...
              kt - boxsize < kt(nv)  );
  isource3=find( abs(kt-kt(nv))==boxsize & ...
              it2 + boxsize > it2(nv)   & ...
              it2 - boxsize < it2(nv)   & ...
              jt + boxsize > jt(nv)   & ...
              jt - boxsize < jt(nv)  );
  isource = [isource1; isource2; isource3];   
  NS = length(isource);

  %% set parameters of the least-squares problem.
  alpha2 = 5;
  obs = [Tobs(nv);Sobs(nv);O18obs(nv);Pobs(nv);Nobs(nv);Oobs(nv)];
  NOBS = length(obs);
     
  %% add an equation for mass conservation.
  obs = [ obs ; 1];

  NS = length(srctmp); % how many sources?

  %% add tapering to the solution.
  % But do it this funny way to fit the matlab lsqnonneg function.
  xbase = (1./NS).*ones(NS,1);
  xbase(end+1) = 0;
  xbase = reshape(xbase,length(xbase),1);
    
  %% Set up the local properties matrix, E.
  Etmp = [T(isource)'; S(isource)'; O18(isource)';P(isource)'; ...
          N(isource)';O(isource)'];
  Etmp(end+1,:) = ones(1,NS);

  %% added next 2 lines for nonconservative tracers.
  Etmp(:,end+1) = zeros(NOBS+1,1);
  Etmp(:,end) = [0 0 0 1 15.5 -170 0];

  obserr =[Terr(nv); Serr(nv); O18err(nv); Perr(nv); ...
                Nerr(nv); Oerr(dsttmp)];
  obserr = [obserr; .001]; % .001 is for mass conservation.
       
  W = diag(1./obserr); % equal to W^(-1/2) in TOCIP, Wunsch 1996.
       
  obsnew = W*obs;
  Enew = W*Etmp;
  if alpha2>0 % then tapering is on.
    Enew = [Enew ; (alpha2./NSlocal).*eye(NSlocal+1)];
    obsnew(end+1:end+length(xbase)) = (alpha2./NSlocal).*xbase;
  end

  x0 = zeros(NS+1,1); % first guess of solution, x ("m" in G&H, 2010)
  [x1,Jtmp,restmp,exflagnew] = lsqnonneg(Enew,obsnew,x0,options); % check to see if Jtmp is reasonable.

  % modify/overwrite the global pathways "A" matrix.
  A(nv,:) = 0;
  A(nv,nv) = -1;
  A(nv, isource) = x1(1:end-1);
  dP = x1(end); % nonconservative source.
end
   
   
