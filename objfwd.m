function c = objfwd(x,A,W,cobs,isfc,iint,noncons,iintdst,iintsrc)
% function c = objfwd(x,A,W,cobs,isfc,iint,noncons,iintdst,iintsrc)
%
% Pass A, W, cobs, isfc, iint, noncons, iintdst, iintsrc,
% Last 2 arguments are optional for the case where pathways are changing.
% Gradient "g" output argument is optional.

Nsfc = length(isfc);
Nint = length(iint);
Nfield = size(A,1);

% use nargin to determine whether paths are changing
if nargin==9
  % if paths are changing, use them to reconstruct A matrix.
  Nm = length(iintdst);
  M = sparse(iintdst,iintsrc,x(end-Nm+1:end),Nfield,Nfield);
  diagsum = sum(M(iint,:)')';
  D = sparse(iint,iint,-diagsum,Nfield,Nfield);
  G = A + M + D;
else
  G = A;  % otherwise pass the A matrix in the arguments.
end

[L, U, P, Q, R] = lu (G) ;
d = zeros(Nfield,1);
d(isfc) = x(1:Nsfc);
if noncons
   d(iint) = x(Nsfc+1:Nsfc+Nint); % add entries to "d" for
                                % nonconservative sources
end

c =  Q * (U \ (L \ (P * (R \ d)))) ; % tracer distribution.

