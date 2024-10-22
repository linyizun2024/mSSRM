function [Q,p] = Qp_gen(m,n)
%% Input variables
% m: number of rows of Q
% n: number of columns of Q

%% Setting of parameters
rho = 0.5;

%% Simulation setup
Sigma_Q = eye(n);
for i = 2:n
    for j = 1:(i-1)
        Sigma_Q(i,j) = rho^abs(i-j);
        Sigma_Q(j,i) = Sigma_Q(i,j);
    end
end
    
%% Generate Q
Q = mvnrnd(zeros(n,1),Sigma_Q,m); % m*n
p = 10*(2*rand(1,n)'-1);

