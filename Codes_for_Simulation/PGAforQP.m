function [x,RE,OFV] = PGAforQP(Q,p,x_init,ITER)
%% Input variables
% Q: Data matrix in the model
% p: Data vector in the model
% Model: 1/2*x'*Q*x-p'*x+iota_{R+}(x)

%% Output variables
% x: solution

%% PGA for Quadratic Programming
x = x_init;
invQp = Q\p;
RE= zeros(ITER,1);
OFV = zeros(ITER,1);
x = max(x-Q\(Q*x)+invQp,0);
beta = 1.99/norm(Q,2);
Qx = Q*x;

for k = 1:ITER
    x_pre = x;
    x = max(x-beta*Qx+beta*p,0);  % PGA
    % x = max(x-Q\Qx+invQp,0);  % PPGA with invHessian preconditioner
    Qx = Q*x;
    RE(k) = norm(x-x_pre,2)/norm(x,2);
    OFV(k) = (1/2)*x'*Qx-p'*x;
end

end
