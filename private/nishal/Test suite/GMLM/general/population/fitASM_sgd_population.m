function [fitASM,fval_log] = fitASM_sgd_population(X,Y,num_su,mask,gamma1,gamma2,lam_start, init_point)


X = gpuArray(X);
Y = gpuArray(Y);

%% Lower dimensional X
X= X(mask,:);

%% 
T =size(X,2);
Nc = size(Y,1);
d = size(X,1);
Ns = num_su;
%% Calculate Sigma
% make sure X is 0 mean ! 

variance = var(X(1,:));
Sigma = diag(variance*ones(d,1));
dt = 1/120;

dt = 1/120;
fval_prev = Inf; fval_prev2=Inf;
epsi=10^-6;
beta = 0.5; % for line search

%% Initialize
if(~isempty(init_point))
    K = init_point.K;
    B = init_point.B;
else
    K = gpuArray(2*(rand(d,Ns)-0.5));
    B = gpuArray(2*(rand(Ns,Nc)-0.5));
end
K_xk = K; B_xk = B;
K_ykp1 = K_xk; B_ykp1 = B_xk;
lambda = lam_start;
%% Update K and B

togo=true;
iter=0;
fval_log=[];


for itime = 1:T
    itime
    eta=0.005;
    
x = X(:,itime);

SU_inp = K'*x;

% find alphas 
alpha = gpuArray(zeros(Ns,1,Nc));
 
 % Find alpha
 fmin = Inf ; 
 Kmin=[];Bmin=[];
 for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B(:,icell),[1,1]));
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000000001;
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alpha(:,:,icell) = alpha_cell;
 end
 
  gradK = 0*K; gradB = 0*B;
 for isu=1:Ns 
 [gradK(:,isu),gradB(isu,:)]=gradientKBj(x,Y(:,itime),B(isu,:),K(:,isu),1,Nc,alpha(isu,:,:),dt); 
 end

 K = K - eta * gradK ; 
 B = B - eta * gradB;
 
 fval_log(itime) = gather(f_val_population(K,B,X,Y,dt));
 
 if(fval_log(itime)<fmin)
 fmin = fval_log;
 Kmin = K;
 Bmin = B;
 end
 
 if(rem(itime,1000)==1)
     close all
        figure;
        plot(fval_log);
        plotSU(K,mask);
        pause(0.1);
 end
 
end
fitASM.K = Kmin;
fitASM.B = Bmin;


figure;
plotSU(K,mask);

end