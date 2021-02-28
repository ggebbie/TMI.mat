function [x]= sq(y)
%function [x]= sq(y)
%
% The "squeeze" matlab command is too long.
% Call it with "sq" instead.
%
% SQUEEZE Remove singleton dimensions.
%    B = SQUEEZE(A) returns an array B with the same elements as
%    A but with all the singleton dimensions removed.  A singleton
%    is a dimension such that size(A,dim)==1.  2-D arrays are
%    unaffected by squeeze so that row vectors remain rows.
% 
%    For example,
%        squeeze(rand(2,1,3))
%    is 2-by-3.
% 
%    See also SHIFTDIM.
%
% G. Gebbie, MIT-WHOI, 26 Jan 2001.

x = squeeze(y);

return
