
function [fitASM,fval_log] = fitASM_EM_Population_sparse_split(X,Y,num_su,mask,gamma1,gamma2,lam_start)

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
K = gpuArray(2*(rand(d,Ns)-0.5));
B = gpuArray(2*(rand(Ns,Nc)-0.5));

K_xk = K; B_xk = B;
K_ykp1 = K_xk; B_ykp1 = B_xk;
lambda_list = lam_start*ones(Ns,1);
%% Update K and B

togo=true;
iter=0;
fval_log=[];
fval_log2=[];
while togo
    iter=iter+1;
 SU_inp = K_ykp1'*X;
 alpha = gpuArray(zeros(Ns,T,Nc));
 
 % Find alpha
 for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B_ykp1(:,icell),[1,T]));
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000000001;
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alpha(:,:,icell) = alpha_cell;
 end
 
 gradK = 0*K; gradB = 0*B;
  K_xk_1=0*K_xk; B_xk_1 = 0*B_xk;
  for isu=1:Ns
 togo_grad=true;
 fval_log=[];
 while togo_grad
 % find gradient of K and B
 [gradK(:,isu),gradB(isu,:)]=gradientKBj(X,Y,B_ykp1(isu,:),K_ykp1(:,isu),T,Nc,alpha(isu,:,:),dt); 
 
  % iterate over lambda
  K_xkp1 = 0*K_ykp1; B_xkp1 = 0*B_ykp1;
  
    lambda = lambda_list(isu);
    togo_linesearch=true;
    fy = f_val_population_su(K_ykp1(:,isu),B_ykp1(isu,:),X,Y,dt,permute(alpha(isu,:,:),[2,3,1]));
    while togo_linesearch
        % proximal L1 step
        K_xkp1(:,isu) = proximalL1(K_ykp1(:,isu) - lambda*gradK(:,isu) , lambda*gamma1);
        B_xkp1(isu,:) = log(proximalL1(exp(B_ykp1(isu,:) - lambda*gradB(isu,:)), lambda*gamma2));
        
        % line search step
        fz = f_val_population_su(K_xkp1(:,isu),B_xkp1(isu,:),X,Y,dt,permute(alpha(isu,:,:),[2,3,1]));
        ftilde_z_ykp1 = ftilde(K_xkp1(:,isu),B_xkp1(isu,:),K_ykp1(:,isu),B_ykp1(isu,:),gradK(:,isu),gradB(isu,:),lambda,fy); % Computations in ftilde screwed up once entries become -Infinity
        % fz, ftilde(z,ykp1)
        if(fz <= ftilde_z_ykp1)
        togo_linesearch = false;
        else
        lambda = beta*lambda;    
        end
    end
    lambda_list(isu) = lambda;
 
    
  % calculate new y_kp1 .. 
    K_xk_1(:,isu) = K_xk(:,isu); B_xk_1(isu,:) = B_xk(isu,:);
    K_xk(:,isu) = K_xkp1(:,isu) ; B_xk(isu,:) = B_xkp1(isu,:);
    
    wk = iter/(iter+3);
    K_ykp1(:,isu) = K_xk(:,isu) + wk *(K_xk(:,isu) - K_xk_1(:,isu));
    B_ykp1(isu,:) = B_xk(isu,:) + wk *(B_xk(isu,:) - B_xk_1(isu,:));
  


fval = fz + gamma1*sum(abs(K_xk(:,isu))) + gamma2*sum(abs(exp(B_xk(isu,:))));
[lambda,fval]
if(abs(fval_prev-fval)<epsi)
    togo_grad=false;
else
    togo_grad=true;
    fval_prev= fval;
end
fval_log = [fval_log;fval];
  end
  
figure;plot(fval_log);title(sprintf('SU: %d',isu));pause(0.05);
  end
 close all
 fval_2 =f_val_population(K_xk,B_xk,X,Y,dt ) +  gamma1*sum(abs(K_xk(:))) + gamma2*sum(abs(exp(B_xk(:))));
 if(abs(fval_prev2 - fval_2)< epsi)
     togo=false;
 else
     togo=true;
     fval_prev2 = fval_2;
 end
 fval_log2 = [fval_log2;fval_2];
 
end

figure;
plotSU(K_xk,mask);
pause(0.1);

fitASM.K = K_xk;
fitASM.B = B_xk;

end