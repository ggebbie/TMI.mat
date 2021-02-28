%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
% load the master file: choose either woce2deg or woce4deg 
% (2 deg == 2 x 2 horiz resolution, 33 vertical levels)
% (4 deg == 4 x 4 horiz resolution, 33 vertical levels)
% (woce == trained on the WOCE Global Hydrographic Climatology,
%          Gouretski and Koltermann 2004)

load A  %% choose between the 4 deg and 2 deg MAT-files
 
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
 Xlon = [12.38 -18.1 -61.25 -14.23 -10.2 118.34 -10.17 51.95 -84];
 % deg E.
 Xlon(Xlon<0) = Xlon(Xlon<0)+360;
 Xlat = [-23.32 12.4 12.08 -44.15 37.8 -9.65 37.8 10.77 -3.6];
 Xdepth = [1967 3223 1299 3770 3146 2100 3166 1580 3210]; % meters.
 
 % Find the coordinate on the grid by linear interpolation. 
 LONtmp = [LON(1)-4 ; LON; LON(end)+4];  %% watch out for the
                                         %wraparound of the grid.
 iX = interp1(LONtmp,0:length(LON)+1,Xlon,'linear')
 jX = interp1(LAT,1:length(LAT),Xlat,'linear')
 kX = interp1(DEPTH,1:length(DEPTH),Xdepth,'linear')
 kX(Xdepth>5500) = 33; % if deeper than the bottom level, put it on
                       % the bottom level.

% $$$  % watch out for wraparound again.
% $$$  i2 = i;
% $$$  i2( i-iX > NX/2) = i2( i-iX > NX/2) - NX;
% $$$  i2( iX-i > NX/2) = i2( iX-i > NX/2) + NX;
  
 % find the closest gridpoints to the location of interest.
 mask = zeros(Nfield,length(iX)); 

 for nc = 1:length(iX)
   [dist,loc] = sort((i-iX(nc)).^2 + (j-jX(nc)).^2 + (k-kX(nc)).^2);
   dist = dist+0.1; % to eliminate singularity if point matches up
                  % exactly with the grid. 
 
   mask(loc(1:6),nc) =(1./dist(1:6))./sum(1./dist(1:6)); % mask adds up
                                                    % to 1.
 
 end
 
 vtot = R' \ (P' * (L' \ (U' \ (Q' * mask)))); 

 G = vtot(k==1,:)';
 y = [-0.43 -0.19 -0.87 0.37 0.43 -0.42 0.64 -0.15 0.61]';
 x = G\y; 
 
 alpha = 0.002;
 x = G'*inv(G*G'+alpha.*eye(9))*y;
 G*x-y
 std(G*x-y)
 std(x)



