function [all_loss,all_v] = PGAforL0QuadProg(Qeps,p,m,initv,ITER)
%% Input variables
% Q: Data matrix in the model
% p: Data vector in the model
% Model: 1/2*v'*Qeps*v-p'*v+iota_m(v)+iota_{R+}(v)

%% PGAforL0QuadProg
v = initv;
n = numel(v);
Qv = Qeps*v;
beta = 0.99/norm(Qeps,2);
all_loss = zeros(ITER,1);
all_v = zeros(n,ITER);

for k = 1:ITER
    v = prox_mpsparse(v-beta*Qv+beta*p,m);
    Qv = Qeps*v;
    all_loss(k) = 1/2*v'*Qv-p'*v;
    all_v(:,k) = v;
end

end
