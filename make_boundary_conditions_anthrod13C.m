%% Make the boundary conditions in the proper form for the 
%  TMI transient tracer simulation.

% Options: 1. Fixed boundary condition (nothing to be done)
%          2. Time varying boundary condition
%
% returns Csfc: the time-evolving boundary conditions field at the surface

switch boundary_type
  case{'anthrod13C'}
    
    % this is the major work for the anthrod13C case. 
    % Need to compute the SST timeseries at each location.
    
    % there are Nsfc surface boundary conditions to fill in.
    Nsfc = length(kt==1);
    for ns = 1:Nsfc
        % 1. find the nearest Devries and Primeau point to get the
        % ETD parameters.
        
        % 1a. Pick the times on which you will get the timeseries.
        tsfc = t0:tf; % in years
        
        % 2. Do a convolution with the ETD and atmospheric d13C.
        
        % 3. Save the d13C_surface timeseries, "Ctimeseries".
        Csfc(ns,:) = Ctimeseries;
    end
  
  case {'fixed'}

    disp('fixed boundary conditions: nothing to be done')


  case {'varying'}

    disp('varying boundary conditions')
    
    Nsfc = sum(kt==1);
    % simple example. ramping up boundary conditions.
    tsfc(1) = 0; %[years]
    tsfc(2) = 50; 
    
    Csfc(1,:) = ones(Nsfc,1); % surface boundary condition at t=tsfc(1)
    Csfc(2,:) = zeros(Nsfc,1);  % surface boundary condition at t=tsfc(2)
end





