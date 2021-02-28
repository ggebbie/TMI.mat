%% Make the boundary conditions in the proper form for the 
%  TMI transient tracer simulation.

% Options: 1. Fixed boundary condition (nothing to be done)
%          2. Time varying boundary condition
%
% returns Csfc: the time-evolving boundary conditions field at the surface

switch boundary_type
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





