%% Make the initial conditions in the proper form for the 
%  TMI transient tracer simulation.

% Examples: 1. Define a surface box to be dyed, 
%           2. Define an arbitrary surface concentration.
% Final step, mix the mixed layer.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Define a surface box.

% define the surface patch by the bounding latitude and longitude.
 lat_lo = 50; % 50 N, for example.
 lat_hi = 60;
 
 lon_lo = -50;
 lon_hi = 0;
 
 loc = find( LAT(j) < lat_hi & LAT(j)> lat_lo & ...
             LON(i) > lon_lo & LON(i)<lon_hi & inmixlyr); 

% Note: potential longitudinal wraparound problems in longitude are not accounted for
% in the "loc" statement.
 
 N = length(i);
 c = zeros(N,1);
 c(loc) = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: DIY some other surface concentration pattern.

% d = some defined surface contribution.

% Now mix the mixed layer.
   c = mixit(d,it,jt,kt,inmixlyr); 
