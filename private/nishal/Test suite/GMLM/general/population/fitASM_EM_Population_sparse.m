
function fitASM = fitASM_EM_Population_sparse(X,Y,num_su,mask,gamma1,gamma2)

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

togo=true;
dt = 1/120;
fval_prev = Inf;
epsi=10^-8;
beta = 0.8; % for line search
while togo
 SU_inp = K'*X;
 alpha = gpuArray(zeros(Nc,Ns,T));
 
 % Find alpha
 for icell=1:Nc
    cell_inp  = exp(SU_inp + repmat(B(:,icell),[1,T]));
    alpha_cell = cell_inp ./ repmat(sum(cell_inp,1),[Ns,1]);
    alpha(icell,:,:) = alpha_cell;
 end
 
 % find gradient of K and B
 gradientKB 
 
  % iterate over lambda
        % proximal L1 step
        
        % line search step
 
  % update K and B
 
    B=Bnew;
    K=Knew;
  
  
  
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


end

figure;
plotSU(K,mask);
pause(0.1);

fitASM.K = K;
fitASM.B = B;

end