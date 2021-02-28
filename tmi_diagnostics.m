%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
%  example 1: track surface waters into interior
%  example 2: Find the volume that has originated from surface
%  example 3: Find surface origin of an interior box
%  example 4: Reconstruct steady-state tracer field.
%  example 5: Construct the pathways matrix from observations.
%  example 6: Follow the adjoint transport pathways.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the master file: choose either woce2deg or woce4deg 
% (2 deg == 2 x 2 horiz resolution, 33 vertical levels)
% (4 deg == 4 x 4 horiz resolution, 33 vertical levels)
% (woce == trained on the WOCE Global Hydrographic Climatology,
%          Gouretski and Koltermann 2004)

load A  %% choose a A*.mat file
 
% variables in A*.mat %%%
%
% A: pathways matrix
% i,j,k: vectors that contain the spatial coordinates 
%        consistent with the A matrix.
% LAT,LON,DEPTH: the latitude, longitude and depth of the 3-d grid. 
% dP = amount of remineralized phosphate in each box.

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
 contourf(LAT,-DEPTH,squeeze(C(:,:,isec)),0:0.05:1) % a sample plot at 22 W.
 
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
  load tracerobs  % choose tracerobs*
  
  dT = Tobs;      % set rhs vector to the observations.
  dT(k>1) = 0;    % no internal sinks or sources.

  % a first guess: observed surface boundary conditions are perfect.  
  Tmod1 =  Q * (U \ (L \ (P * (R \ dT))));
  Tmod1_field = vector_to_field(Tmod1,i,j,k);
 
  % get cost function (J) for first guess temperature field.
  % here the data-model misfit is weighted by the expected error. 
  Nfield = length(Tobs);

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

    tmp1 = invWT*H2;
    HT = H2'*tmp1;
    fT = 2.*H1'*tmp1;
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
% Example 5: How to modify or reconstruct the TMI pathways,            %
%            i.e., the A matrix.                                       %
%                                                                      %
% Each row of A is determined by the solution of a non-negative        %
% least squares problem via the method of Lawson and Hanson.           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
it = i; jt = j; kt = k;
NX = max(it); NY = max(jt); NZ = max(kt); Nfield = size(A,1);

% save truth to compare later.
Atruth = A;
dPtruth = dP;

for nv = 1:Nfield

  % If in mixed layer, then stop. The local inversion is not applicable.
  if ~inmixlyr(nv)
  
    nv
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
    alpha2 = 0;
    obs = [Tobs(nv);Sobs(nv);O18obs(nv);Pobs(nv);Nobs(nv);Oobs(nv)];
    NOBS = length(obs);
     
    %% add an equation for mass conservation.
    obs = [ obs ; 1];

    %% add tapering to the solution.
    % But do it this funny way to fit the matlab lsqnonneg function.
    xbase = (1./NS).*ones(NS,1);
    xbase(end+1) = 0;
    xbase = reshape(xbase,length(xbase),1);
    
    %% Set up the local properties matrix, E.
    Etmp = [Tobs(isource)'; Sobs(isource)'; O18obs(isource)';Pobs(isource)'; ...
            Nobs(isource)';Oobs(isource)'];
    Etmp(end+1,:) = ones(1,NS);

    %% added next 2 lines for nonconservative tracers.
    Etmp(:,end+1) = zeros(NOBS+1,1);
    Etmp(:,end) = [0 0 0 1 15.5 -170 0];

    obserr =[Terr(nv); Serr(nv); O18err(nv); Perr(nv); ...
                Nerr(nv); Oerr(nv)];
    obserr = [obserr; .001]; % .001 is for mass conservation.
       
    W = diag(1./obserr); % equal to W^(-1/2) in TOCIP, Wunsch 1996.
       
    obsnew = W*obs;
    Enew = W*Etmp;
    if alpha2>0 % then tapering is on.
      Enew = [Enew ; (alpha2./NS).*eye(NS+1)];
      obsnew(end+1:end+length(xbase)) = (alpha2./NS).*xbase;
    end

    x0 = zeros(NS+1,1); % first guess of solution, x ("m" in G&H, 2010)
    [x1,Jtmp,restmp,exflagnew] = lsqnonneg(Enew,obsnew); % check to see if Jtmp is reasonable.

    % modify/overwrite the global pathways "A" matrix.
    A(nv,:) = 0;
    A(nv,nv) = -1;
    A(nv, isource) = x1(1:end-1)./sum(x1(1:end-1)); % normalize just
                                                   % in case.
    dP(nv) = x1(end); % nonconservative source.
  end
end

% compare A to Atruth if using TMI tracer observations.
% compare dP to dPtruth. 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6
%
% Track the waters from an interior point upstream to the surface.
% For each interior location, determine the fraction of waters 
% from location "pt" that go through each interior location before
% getting to the surface.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfield = size(A,1);
d = zeros(Nfield,1);
%pt = maxloc(9); an example southern ocean point from
%watermass-unmixing study.
d(pt) = 1;
A_orig = A;
vsave = nan(Nfield,1);
%%% do a depth profile. 
%% IT SEEMS TO WORK.
%matlabpool open 8
parfor nv = 1:Nfield
  %tic
  nv
    
  A = A_orig; % reset A.
  A(nv,:) = 0;
  A(nv,nv) = 1;
  [L, U, P, Q, R] = lu (A);
  vtot = R' \ (P' * (L' \ (U' \ (Q' * d))));  
  vsave(nv) = vtot(nv);
  %toc
end

Vsave = vector_to_field(vsave,it,jt,kt);
%contourf(LON,LAT,Vtot)

contourf(LAT,-DEPTH,log10(sq(Vsave(:,:,78))))
[c,h]=contourf(LAT,-DEPTH,log10(sq(Vsave(:,:,78))),-8:0.25:0)
[c,h]=contourf(LON,LAT,log10(sq(Vsave(23,:,:))),-8:.25:0)


