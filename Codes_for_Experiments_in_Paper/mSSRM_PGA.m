function [w] = mSSRM_PGA(Param,matR,vecmu)
[T,N] = size(matR);
RE=100;
eI=Param.eps*eye(N);
p=vecmu;
Q=(1/sqrt(T-1))*(matR-(1/T)*ones(T,T)*matR);
QeI=Q'*Q+eI;
%% iteration
Param.alpha=0.999/norm(QeI,2);
w=vecmu;
k = 1;
while k<=Param.iternum && RE(k)>Param.tol
    w1=w;
    w_pre=w-Param.alpha*(QeI*w-p) ;
    w_pre(w_pre<0)=zeros(sum(w_pre<0),1);
    [~,itw]=sort(w_pre,'descend');
    w=zeros(N,1);
    w(itw(1:Param.m))=w_pre(itw(1:Param.m));
    RE(k+1) = norm(w-w1,2)/norm(w1,2);
    k = k+1;
end

if sum(w)==0
    w = zeros(N,1);
else
    w=w/sum(w);
end

end
