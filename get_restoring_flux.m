function [f] = get_restoring_flux(t,C,Cb,tau)
% function [f] = get_restoring_flux(t,C,Cb,tau)
%
% f = forcing (units C/yr)
%
% C = C(t) current tracer value
% Cb(1) target value at t=0
% Cb(2) target value at t=1
% tau = restoring timescale (yr^-1)

f = -(C-Cb)./tau; 










