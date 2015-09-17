
function [fitASM,fval_log] = fitASM_EM_Population_sparse_linearb(X,Y,num_su,mask,gamma1,gamma2,lam_start)

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
fval_prev = Inf;
epsi=10^-8;
beta = 0.5; % for line search

%% Initialize
K = gpuArray(2*(rand(d,Ns)-0.5));
B = exp(gpuArray(2*(rand(Ns,Nc)-0.5)));

K_xk = K; B_xk = B;
K_ykp1 = K_xk; B_ykp1 = B_xk;
lambda = lam_start;
%% Update K and B

togo=true;
iter=0;
fval_log=[];
while togo
    iter=iter+1;
    
    
    if(rem(iter,50)==1)
 SU_inp = exp(K_ykp1'*X);
 
 
 % Find alpha k
 alphak = gpuArray(zeros(Ns,T,Nc));
 for icell=1:Nc
    cell_inp  = (SU_inp).* repmat(B_ykp1(:,icell),[1,T]);
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000000001;
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alphak(:,:,icell) = alpha_cell;
 end
 
 % Find alpha b
 alphab = gpuArray(zeros(Ns,T,Nc));
 for icell=1:Nc
    cell_inp  = (SU_inp).* repmat(B_ykp1(:,icell),[1,T]);
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000000001;
    alpha_cell = (SU_inp) ./ repmat(sum(cell_inp,1),[Ns,1]);
    alphab(:,:,icell) = alpha_cell;
 end
    end
 
 
 % find gradient of K and B
 gradK = 0*K; gradB = 0*B;
 for isu=1:Ns
 [gradK(:,isu),gradB(isu,:)]=gradientKBj_linearb(X,Y,B_ykp1(isu,:),K_ykp1(:,isu),T,Nc,permute(alphak(isu,:,:),[2,3,1]),permute(alphab(isu,:,:),[2,3,1]),dt); 
 end
 
  % iterate over lambda
    togo_linesearch=true;
    fy = f_val_population_linearb(K_ykp1,B_ykp1,X,Y,dt);
    while togo_linesearch
        % proximal L1 step
        K_xkp1 = proximalL1(K_ykp1 - lambda*gradK , lambda*gamma1);
        B_xkp1 = (proximalL1((B_ykp1 - lambda*gradB), lambda*gamma2));
        
        % line search step
        fz = f_val_population_linearb(K_xkp1,B_xkp1,X,Y,dt)
        ftilde_z_ykp1 = ftilde(K_xkp1,B_xkp1,K_ykp1,B_ykp1,gradK,gradB,lambda,fy); % Computations in ftilde screwed up once entries become -Infinity
        % fz, ftilde(z,ykp1)
        if(fz <= ftilde_z_ykp1)
        togo_linesearch = false;
        else
        lambda = beta*lambda;    
        end
    end
    
    
  % calculate new y_kp1 .. 
    K_xk_1 = K_xk; B_xk_1 = B_xk;
    K_xk = K_xkp1 ; B_xk = B_xkp1;
    
    wk = iter/(iter+3);
    K_ykp1 = K_xk + wk *(K_xk - K_xk_1);
    B_ykp1 = B_xk + wk *(B_xk - B_xk_1);
  


fval = fz + gamma1*sum(abs(K_xk(:))) + gamma2*sum(abs((B_xk(:))));
[lambda,fval]
if(abs(fval_prev-fval)<epsi)
    togo=false;
else
    togo=true;
    fval_prev= fval;
end
fval_log = [fval_log;fval];

% if(rem(iter,40)==0)
%  keyboard
% end

end

figure;
plotSU(K_xk,mask);
pause(0.1);

fitASM.K = K_xk;
fitASM.B = B_xk;

end