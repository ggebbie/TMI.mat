function F = get_F(A,Gamma,E)

% FDIC
if size(Gamma,2) <= size(E,1)
  tic;
  tmp0 = A\Gamma;
  F =  E*tmp0;
  clear tmp0
  toc;
else
  tic;
  tmp0 = A'\E';
  F = transpose(Gamma'*tmp0);
  toc;
end
