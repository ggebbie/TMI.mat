function H = get_hessian(G,A,E,W,S);

if size(G,2) <= size(E,1)
  tic;
  tmp0 = A\G;
  EAG =  E*tmp0;
  clear tmp0
  toc;
else
  tic;
  tmp0 = A'\E';
  EAG = transpose(G'*tmp0);
  toc;
end

tmp1 = W*EAG;
tmp2 = EAG'*tmp1;

if isscalar(S)
  H    = tmp2 + S.*speye(size(tmp2,1));
else
  H    = tmp2 + S;
end
