%% Make the initial conditions in the proper form for the 
%  TMI transient tracer simulation.

% Options: 1. Choose a pre-defined surface region.
%          2. Define a surface rectangle to be dyed. 
%
% returns c0: the initial conditions field

switch initial_type
  case {'predefined'}

    disp('predefined-region initial conditions')

    % choose one of my pre-defined surface regions:
    % 1) GLOBAL,2) ANT, 3 SUBANT, 4 NATL, 5 NPAC, 
    % 6) TROP, 7. ARC, 8 MED, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
    % 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
    % 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP

    wm = 1; % choose 1 to 19

    % The next lines use the pre-defined regions. 
    load c_all_4deg
    c0 = squeeze(c_all(:,wm));

  case {'rectangle'}

    disp('rectangle initial conditions')
    %% Define a surface rectangle.

    % define the surface patch by the bounding latitude and longitude.
    lat_lo = 50; % 50 N, for example.
    lat_hi = 60;
 
    lon_lo = -50;
    lon_hi = 0;

    %% take care of the longitudinal wraparound.
    if lon_lo<0
      lon_lo = lon_lo + 360;
    end
    if lon_hi<0
      lon_hi = lon_hi + 360;
    end
    
    if lon_lo > lon_hi % wraparound the prime meridian
       loc = find( LAT(jt) <= lat_hi & LAT(jt) >= lat_lo & ...
                  (LON(it) >= lon_lo | LON(it) <= lon_hi) & inmixlyr); 
    else
       loc = find( LAT(jt) <= lat_hi & LAT(jt)>= lat_lo & ...
                   LON(it) >= lon_lo & LON(it)<= lon_hi & inmixlyr); 
    end

    N = length(it);
    c0 = zeros(N,1);
    c0(loc) = 1;

end
 

% Now mix the mixed layer.
%c = mixit(d,it,jt,kt,inmixlyr); 




