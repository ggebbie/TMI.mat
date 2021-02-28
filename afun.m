function y = afun(x,Gamma,A,E,W,S);

ETWE = E'*(W*E);

tmp1 = A\(Gamma*x);
tmp2 = ETWE*tmp1;
tmp3 = A'\tmp2;

y    = Gamma'*tmp3 + S*x;
