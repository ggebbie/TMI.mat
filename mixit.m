function [c] = mixit(d,it,jt,kt,inmixlyr);
% function [c] = mixit(d,it,jt,kt,inmixlyr);
%
% Take the distribution, d, defined only at the surface,
% and return the distribution, c, that has a mixed layer. 

Nfield = length(d);
Nsfc = sum(kt==1);
isfc = find(kt==1);

c = 0.*d;

for ns = 1:Nsfc
  c( it == it(isfc(ns)) & jt==jt(isfc(ns)) & inmixlyr) = d(isfc(ns));
end

  
