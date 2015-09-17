
function fitASM = fitASM_EM_Population(X,Y,num_su,mask)

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


%% Initialize
K =  2*(rand(d,Ns)-0.5);
B = 2*(rand(Ns,Nc)-0.5);




%% Update K and B

togo=true;
dt = 1/120;
fval_prev = Inf;
epsi=10^-8;
while togo
 SU_inp = K'*X;
 alpha = gpuArray(zeros(Nc,Ns,T));
 
 % Find alpha
 for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B(:,icell),[1,T]));
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alpha(icell,:,:) = alpha_cell;
 end
 
 % find Knew
 Knew = 0*K;
 fval=0;
    for isu=1:Ns
        [fval_su,Knew(:,isu)] = update_Kj_ASM_population(X,Y,B(isu,:),K(:,isu),T,Nc,SU_inp,squeeze(alpha(:,isu,:)),dt);
    fval = fval + fval_su;
    end
    
 % find Bnew   
 Bnew = 0*B;
 total_SU_act = sum(exp(Knew'*X),2);
 for isu=1:Ns
 Bnew(isu,:) = Bnew(isu,:) + gather(sum(squeeze(alpha(:,isu,:))'.*Y,2)'/T);
 end
Bnew =log(Bnew);
Bnew = Bnew - log(repmat(total_SU_act,[1,Nc]));
Bnew = Bnew - log(dt/T);

% calulate negative LL 
SU_inp_new = Knew'*X;
LL=0;
for icell = 1:Nc
    cell_inp = exp(SU_inp_new + repmat(Bnew(:,icell),[1,T]));
    LL = LL + sum(cell_inp(:))*dt/T - (1/T)*log(sum(cell_inp,1))*Y(icell,:)';
end
fval = LL


if(abs(fval_prev-fval)<epsi)
    togo=false;
else
    togo=true;
    fval_prev= fval;
end

B=Bnew;
K=Knew;

end

figure;
plotSU(K,mask);
pause(0.1);

fitASM.K = K;
fitASM.B = B;

end