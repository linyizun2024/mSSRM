function [CW,sharpe] = run_mSSRM_PGA(win_size,data,m)
Param.winsize=win_size;
Param.m=m;
Param.iternum=1e4;
Param.tol=1e-5;
Param.eps=1e-3;

fullR = (data-1);
[fullT,N] = size(fullR);
T_end = fullT;
all_w = ones(N,fullT)/N;
CW = zeros(T_end,1);
S = 1;

for t = 1:T_end
    if t>5
        if t<=Param.winsize
            win_start = 1;
        else
            win_start = t-Param.winsize;
        end
        win_end = t-1;
        T = win_end-win_start+1;
        matR = fullR(win_start:win_end,:);
        vecmu = (sum(matR)/T)';
        [w] = mSSRM_PGA(Param,matR,vecmu);
        all_w(:,t) = w;
        if sum(w)~=0
            S = S*data(t,:)*all_w(:,t);
        end
    end   
    CW(t) = S;
    if mod(t,10)==0
        fprintf('***This is the %dth trade day***\n\n',t);
    end
end

A=tick2ret(CW);
a=mean(A);
b=std(A);
sharpe=a/b;

end