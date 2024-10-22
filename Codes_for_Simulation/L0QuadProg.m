function [opt] = L0QuadProg(Qeps,p,m,ITER)
%% Input variables
% Qeps: Data matrix in the model
% p: Data vector in the model
% m: sparsity
% Model: 1/2*v'*Q*v-p'*v+iota_m(v)+iota_{R+}(v)

%% Output variables
% v: global minimizer

%% Brute-force method for L0 regression
opt.loss = inf;
n = size(Qeps,2);
opt.supp = zeros(1,n);
opt.v = zeros(n,1);
allcomb = combnk(1:n,m);
x_init = zeros(m,1);
for i = 1:size(allcomb,1)
    support = allcomb(i,:);
    Qeps_supp = Qeps(support,support);
    p_supp = p(support);
    [x,RE_PGA,OFV_PGA] = PGAforQP(Qeps_supp,p_supp,x_init,ITER);
    v = zeros(n,1);
    v(support) = x;
    lossvalue = (1/2)*v'*Qeps*v-p'*v;
    if lossvalue<opt.loss
        opt.loss = lossvalue;
        opt.supp = support;
        opt.v = v;
    end
end

end
