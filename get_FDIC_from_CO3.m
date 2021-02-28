function [FDIC] = get_FDIC_from_CO3(A,Gamma)

% get EDIC


% FDIC
if size(Gamma,2) <= size(EDIC,1)
  tic;
  tmp0 = A\Gamma;
  FDIC =  EDIC*tmp0;
  clear tmp0
  toc;
else
  tic;
  tmp0 = A'\EDIC';
  FDIC = transpose(Gamma'*tmp0);
  toc;
end
