clear;
close all;
clc;
%% Setting of parameters
testnum = 10000;
epsilon = 0.001;
M = 50;
N = 10;
m = 3;
initv = ones(N,1);%/N;
%initv = zeros(N,1);
ITER_Opt = 100;
ITER_PGA = 500;
finalEFV = zeros(1,testnum);
finalEIS = zeros(1,testnum);

for k = 1:testnum
    %% Generate the data for regression
    [Q,p] = Qp_gen(M,N);
    Qeps = Q'*Q+epsilon*eye(size(Q,2));
    
    %% Call regression functions
    [opt0] = L0QuadProg(Qeps,p,m,ITER_Opt);
    [all_loss,all_v] = PGAforL0QuadProg(Qeps,p,m,initv,ITER_PGA);
    
    NEFV = abs(all_loss-opt0.loss)/abs(opt0.loss);  % Normalized Error of Function Value;
    NEIS = zeros(ITER_PGA,1);  % Normalized Error of Iterative Sequence;
    normoptv = norm(opt0.v,2);
    for i = 1:ITER_PGA
        NEIS(i) = (norm(all_v(:,i)-opt0.v,2))/normoptv;
    end
    finalEIS(k) = NEIS(ITER_PGA);
    finalEFV(k) = NEFV(ITER_PGA);
end

globnum1 = sum(finalEIS<1E-10);
globnum2 = sum(finalEFV<1E-10);

nonglobnum1 = sum(finalEIS>=1E-5);
nonglobnum2 = sum(finalEFV>=1E-5);

avererror1 = sum(finalEIS(finalEIS>=1E-5))/nonglobnum1;
avererror2 = sum(finalEFV(finalEFV>=1E-5))/nonglobnum2;

maxNEIS = max(finalEIS);
maxNEFV = max(finalEFV);

% data = [finalEIS;finalEFV];
% writematrix(data, 'output.xlsx');

figure;
plot(1:ITER_PGA,NEFV,'-r','LineWidth',3);
xlabel('Iteration number');ylabel('$|f(v^k)-f(v^*)|/|f(v^*)|$','Interpreter','latex');
set(gca,'FontSize',30)
% 
figure;
plot(1:ITER_PGA,NEIS,'-b','LineWidth',2);
xlabel('Iteration number');ylabel('$\|v^k-v^*\|_2/\|v^*\|_2$','Interpreter','latex');
set(gca,'FontSize',30)

% savefig = 0;
% if savefig == 1
%     set(gcf,'Units','Inches');
%     pos = get(gcf,'Position');
%     set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
%     print(gcf,'myplot','-dpdf','-r0')
% end