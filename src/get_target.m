function [Cb] = get_target(t,tsfc,Csfc)
% function [Cb] = get_target(t,tsfc,Csfc)
%
% Cb = target surface field
%
% t = current time
% tsfc = surface boundary input time.
% Csfc = surface boundary input field.


% linear interp.
tinterp = interp1(tsfc,1:length(tsfc),t);

Cb = (ceil(tinterp)-tinterp).*Csfc(floor(tinterp),:) + (tinterp-floor(tinterp)).*Csfc(ceil(tinterp),:);
