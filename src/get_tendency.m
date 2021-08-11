function [dCdt] = get_tendency(t,C,L,B,isfc,tsfc,Csfc,tau)
% function [dCdt] = get_tendency(t,C,L,B,isfc,tsfc,Csfc,tau)

display('get tendency: t=')
t
dCdt = L*C; 

if nargin == 8
  Cb = get_target(t,tsfc,Csfc);
  f = get_restoring_flux(t,C(isfc),Cb',tau);
  dCdt = dCdt + B*f;
end









