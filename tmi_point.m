function tmi_point(lon,lat,depth,A)
% tmi_point
%
% automate all diagnostics at a point (lon,lat,depth)
% A is the pathways matrix.

 [L, U, P, Q, R] = lu (A) ;

 % Also, can read in predefined surface patches in the file
 % d_all.mat, where the surface patches are defined as
 % oceanographically-relevant regions: 
 % 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
 % 6) ARC, 7 MED, 8) TROP, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
 % 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
 % 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP
 %
 % e.g.,
 load d_all_4deg; 
 d = d_all_4deg;
 
 % to get quantity of dyed water throughout ocean:
 c = Q * (U \ (L \ (P * (R \ d))));  
 C = vector_to_field(c,i,j,k);  % g is the 3-d field of water-mass concentration.

 %% figure: profile of major water masses through the point of
 %  interest.
 
 
 lon_section = 330;
 isec = find(LON==330);
 contourf(LAT,-DEPTH,squeeze(C(:,:,isec)),0:0.05:1) % a sample plot at 22 W.

 
 NZ = max(k);
 NY = max(j);
 NX = max(i);
 Nfield = size(A,1);
 
 
 % Find the coordinate on the grid by linear interpolation. 
 LONtmp = [LON(1)-4 ; LON; LON(end)+4];  %% watch out for the
                                         %wraparound of the grid.
 iX = interp1(LONtmp,0:length(LON)+1,lon,'linear')
 jX = interp1(LAT,1:length(LAT),lat,'linear')
 kX = interp1(DEPTH,1:length(DEPTH),depth,'linear')
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

 figure
 % contourf(LON,LAT,Vtot,[0 .001 .005 .01 .05 .1 .5]) 
 contourf(LON,LAT,log10(Vtot)) 

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6
%
% Track the waters from an interior point upstream to the surface.
% For each interior location, determine the fraction of waters 
% from location "pt" that go through each interior location before
% getting to the surface.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 2
if method == 1
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
elseif method==2
  AT = A';
  AT(
  
  
end

Vsave = vector_to_field(vsave,it,jt,kt);
%contourf(LON,LAT,Vtot)

contourf(LAT,-DEPTH,log10(sq(Vsave(:,:,78))))
[c,h]=contourf(LAT,-DEPTH,log10(sq(Vsave(:,:,78))),-8:0.25:0)
[c,h]=contourf(LON,LAT,log10(sq(Vsave(23,:,:))),-8:.25:0)


