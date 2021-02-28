function [J,g] = objfun(x,A,W,cobs,isfc,iint,noncons,iintdst,iintsrc)
% function [J,g] = objfun(x,A,W,cobs,isfc,iint,noncons,iintdst,iintsrc)
%
% Pass A, W, cobs, isfc, iint, noncons, iintdst, iintsrc,
% Last 2 arguments are optional for the case where pathways are changing.
% Gradient "g" output argument is optional.
nargin
adjoint = (nargout == 2) % calculate gradient with adjoint.
mod_paths = (nargin==9) % interior pathways = control
mod_noncons = max(abs(noncons))>0 % interior nonconservative
                                   % source = control

Nx = length(x);
Nsfc = length(isfc);
Nint = length(iint);
Nfield = size(A,1);
Nc = length(noncons);

J = 0;
g = zeros(Nx,1);

if mod_paths
  % if paths are changing, use them to reconstruct A matrix.
  Nm = length(iintdst);
  inmixlyr = setxor(1:Nfield,iint); % used in adjoint.
  Y = sparse(iintdst,1:length(iintsrc),ones(Nm,1),Nfield,Nm);
  Y = Y(iint,:);
  im = length(x)-Nm+1:length(x);
  m = x(im);
  diagsum = Y*m;
  D = sparse(iint,iint,-diagsum,Nfield,Nfield);
  M = sparse(iintdst,iintsrc,m,Nfield,Nfield);
  G = A + M + D;

  J = J + 1./(5/6).^2. * (m-(1/6))'*(m-(1/6)); % expect m= 1/6 +- 5/6 in absence of
                                  % other info.
  J = J + 1./(0.1).^2.*(Y*m-ones(Nint,1))'*(Y*m-ones(Nint,1));     % sum(m) has
                                               % uncertainty
                                               % 0.1. Purely
                                               
  % non-negative barrier function
  mlo = find(m<1e-4);
  J = J + 1e8.*sum((m(mlo)-1e-4).^4);

  if adjoint
    g(im) =  2./(5/6).^2 .* (m-(1/6)); % expect m= 1/6 +- 5/6 in absence of
                                  % other info.
    g(im) = g(im) + 2./(0.1).^2*Y'*(Y*m-ones(Nint,1));     % sum(m) has uncertainty 0.1. Purely
                             % numerical device.

    g(im(mlo)) = g(im(mlo)) + 4e8.*(m(mlo)-1e-3).^3;
  end
else
  G = A;  % otherwise pass the A matrix in the arguments.
end

[L, U, P, Q, R] = lu (G) ;

if mod_noncons
  iq = Nc*Nsfc+1:Nc*Nsfc+Nint;
  q = x(iq);
  J = J + 1./(0.3.^2)*q'*q;
  
  qlo = find(q<1e-4);
  J = J + 1e8.*sum((q(qlo)-1e-4).^4);

  if adjoint
    g(iq) = 2./(0.3^2) * q; % q is expected to be 0.3 or smaller.
                         % a calculation for the adjoint later.
    g(iq(qlo)) = g(iq(qlo)) + 4e8.*(q(qlo)-1e-4).^3;
  end
end

% now scroll through the tracers.
for nc = 1:Nc
  nc
  d = zeros(Nfield,1);
  isfctmp = ((nc-1)*Nsfc)+1:nc*Nsfc;
  d(isfc) = x(isfctmp);
  if noncons(nc)~=0 
    d(iint) = -noncons(nc).*q; % add entries to "d" for
                                % nonconservative sources
  end

  c =  Q * (U \ (L \ (P * (R \ d)))) ; % tracer distribution.

  misfit = c-cobs(:,nc);
  J = J+misfit'*W{nc}*misfit
  (misfit'*W{nc}*misfit)./Nfield

  %% get Gradient of cost function w.r.t. independent parameters.
  if nargout == 2
    dc = 2*W{nc}*misfit; % forcing by tracer misfit.
  
    gc = R' \ (P' * (L' \ (U' \ (Q' * dc)))); 
    g(isfctmp) = g(isfctmp) + gc(isfc);
  
    if noncons(nc)~=0
      g(iq) = g(iq)-gc(iint).*noncons(nc); % plus add normal sensitivity term.
    end
  
    if mod_paths % get sensitivity to pathways
      csrc = c(iintsrc)-c(iintdst);
      dAdm = sparse(1:length(iintsrc),iintdst,csrc,Nm,Nfield); 
      dAdm(:,inmixlyr) = [];
    
      g(im) = g(im)-dAdm*gc(iint);
    end % pathways loop
  end % adjoint loop
end % tracers loop

   
