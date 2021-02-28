function [fval, exitflag, output, x] = hessianprob(A,G,W,y,u0,d,l,u);
 % tmi  
 % Get xstart, u, l, B, A, f
   
mtxmpy = @get_hessian; % function handle to qpbox4mult nested subfunction

% Choose the HessMult option
% Override the TolPCG option
options = optimset('HessMult',mtxmpy,'TolPcg',0.01,'Display','iter','MaxIter',10,'MaxPCGIter',10);

lb =  l.*ones(length(u0),1)-G'*d;
ub =  u.*ones(length(u0),1)-G'*d;

[L, U, P, Q, R] = lu (A) ;
Nfield = size(A,1);
r = Q * (U \ (L \ (P * (R \ d)))) - y;  
f =  W*r;

f = R' \ (P' * (L' \ (U' \ (Q' * f))));  % inverse of transpose

%clear L U P Q R
f = G'*f;

I = speye(length(u0));
[x, fval, exitflag, output]=quadprog(I,f,[],[],[],[],lb,ub,u0,options);

    function H = get_hessian(I,x); % no conditioning, i.e. "I"
       xsave = x(:,1);
       x = G*x;
       x = Q * (U \ (L \ (P * (R \ x))));  % inverse multiply
       x = W*x;
       x = R' \ (P' * (L' \ (U' \ (Q' * x))));  % inverse of
                                                % transpose then multiply

       H = G'*x;

        J = (xsave'*H + f'*xsave + r'*(W*r))
        xsave'*H 
        f'*xsave
        r'*(W*r)
    end
end
    
