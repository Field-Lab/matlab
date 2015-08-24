
function fitASM = fitASM_MEL_EM_Population(X,Y,num_su,mask)

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
K = 0.01*randn(d,Ns);
B = 0.01*randn(Ns,Nc);




%% Update K and B

LL_prev = Inf;
eps = 10^-8;
togo=true;
icnt=0;
K_log = cell(100,1);
B_log = cell(100,1);
LLmin=Inf; imin=0; max_iter=20;

while togo
    icnt=icnt+1;
SU_inp = K'*X;
Ya = zeros(Ns,1);
K_sum = gpuArray(zeros(d,Ns));
Bnew = zeros(Ns,Nc);
Ya_log = gpuArray(zeros(Ns,Nc));
LL=0; LL_expec=0;
for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B(:,icell),[1,T]));
    alpha = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    Ya_log(:,icell) = alpha*(Y(icell,:)'/T);
    Ya = Ya + Ya_log(:,icell); 
  
    
    for isu=1:Ns
        K_sum(:,isu) = K_sum(:,isu) + X*(alpha(isu,:).*Y(icell,:)/T)';
    end
    
    LL = LL + (1/T)*(sum(cell_inp(:))*dt  - log(sum(cell_inp,1))*Y(icell,:)'); 
    LL_expec =  LL_expec + (sum(exp(0.5*diag(K'*Sigma*K)).*exp(B(:,icell)))*dt - log(sum(cell_inp,1))*Y(icell,:)')/T;
end

Knew =Sigma\( K_sum./repmat(Ya',[d,1]));

KSigK =  -0.5*diag(Knew'*(Sigma * Knew)); 
Bnew = repmat(KSigK,[1,Nc])   +log(Ya_log) -log(dt);

if(abs(LL - LL_prev)<eps)
    togo=false;
    
else
    togo=true & icnt<max_iter;
    LL_prev=LL
    LL_expec
    
end

    K=Knew;
    B=Bnew;
%     figure;
%     hist(K(:));
%     pause(0.1);
plotSU(K,mask);

K_log{icnt} = gather(K);
B_log{icnt}=gather(B);

if(LL<LLmin)
    imin = icnt;
    LLmin=LL;
end 

close all


end

fitASM.K = K_log{imin};
fitASM.B = B_log{imin};

plotSU(K_log{imin},mask);pause(0.1);




end