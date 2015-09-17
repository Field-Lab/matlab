% fit ASM population - ADMM 
function [fitASM,fval_log2] = fitASM_EM_Population_sparse_split_admm(X,Y,num_su,mask,gamma1,gamma2,lam_start,initVal)

X = gpuArray(X);
Y = gpuArray(Y);

%% Lower dimensional X
X= X(mask,:);

%% 
T =size(X,2);
Nc = size(Y,1);
d = size(X,1);
Ns = num_su;
scale_cell=sum(Y,2); % experiment with this!
%% Calculate Sigma
% make sure X is 0 mean ! 

variance = var(X(1,:));
Sigma = diag(variance*ones(d,1));
dt = 1/120;

dt = 1/120;

epsi=10^-6;

rho =lam_start; % for ADMM 

%% Initialize

if(isempty(initVal))
K = gpuArray(2*(rand(d,Ns)-0.5));
B = gpuArray(2*(rand(Ns,Nc)-0.5));
else
K = initVal.K;
B = initVal.B;
end

fval_prev2=f_val_population(K,B,X,Y,dt,scale_cell) +  gamma1*sum(abs(K(:))) + gamma2*sum(abs(exp(B(:))));
K_xk = K; B_xk = B;


%% Update K and B

togo=true;
iter=0;
fval_log=[];
fval_log2=[fval_prev2];
while togo
    iter=iter+1;
 SU_inp = K'*X;
 alpha = gpuArray(zeros(Ns,T,Nc));
 
 % Find alpha
 for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B(:,icell),[1,T]));
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000000001;
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alpha(:,:,icell) = alpha_cell;
 end
 
 for isu=1:Ns
 
 togo_grad=true;
 fval_log=[];
 K_zk = K(:,isu); B_zk=B(isu,:);
 u_Kk = 0*K_zk; u_Bk = 0*B_zk;
 fval_prev = Inf;
 while togo_grad
 % minimize convex upper bound
 
 [K_xkp1,B_xkp1] = update_Kj_ASM_population_admm (K_zk,B_zk,u_Kk,u_Bk,rho,X,Y,permute(alpha(isu,:,:),[2,3,1]),dt,scale_cell);
 
 % Project onto L1 set
 % [gradK(:,isu),gradB(isu,:)]=gradientKBj(X,Y,B_ykp1(isu,:),K_ykp1(:,isu),T,Nc,alpha(isu,:,:),dt); 
 
        K_zkp1 = proximalL1(K_xkp1 + u_Kk , gamma1/rho);
        %B_zkp1 = log(proximalL1(exp(B_xkp1) - u_Bk, gamma2/rho));
        B_zkp1 = B_xkp1;
        
        
 % update u
        u_Kk = u_Kk + (K_xkp1 - K_zkp1);
        u_Bk = u_Bk + (B_xkp1 - B_zkp1);
        
        K_zk = K_zkp1; B_zk = B_zkp1;


fval = norm((K_xkp1 - K_zkp1)) + norm(B_xkp1 - B_zkp1);%f_val_population_su(K_xkp1,B_xkp1,X,Y,dt,permute(alpha(isu,:,:),[2,3,1])) + gamma1*sum(abs(K_xk(:,isu))) + gamma2*sum(abs(exp(B_xk(isu,:))));
if(abs(fval_prev-fval)<epsi)
    togo_grad=false;
    K(:,isu) = K_xkp1; B(isu,:) = B_xkp1;
else
    togo_grad=true;
    fval_prev= fval;
end
fval_log = [fval_log;fval];
  end
  
%figure;plot(log(fval_log));title(sprintf('SU: %d',isu));pause(0.05);
 end
  
 close all
 fval_2 =f_val_population(K,B,X,Y,dt,scale_cell) +  gamma1*sum(abs(K(:))) + gamma2*sum(abs(exp(B(:))))
 if(abs(fval_prev2 - fval_2)< epsi)
     togo=false;
 else
     togo=true;
     fval_prev2 = fval_2;
 end
 fval_log2 = [fval_log2;fval_2];
 
end

%figure;
%plotSU(K,mask);
%pause(0.1);

fitASM.K = K;
fitASM.B = B;

end