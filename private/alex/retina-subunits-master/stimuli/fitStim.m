function prsEst = fitStim(fit1,fit2,testMode)

% set up an fmincon routine
initPrs = [normMat(vector(randn(1,fit1.C)),'L2',1)];

% create handle to objective function
fun = @(prs) fitStim_errFun(prs(:),fit1,fit2);
x0 = initPrs(:);

if testMode == 1
    DerivCheck_Elts(fun,x0);
    return
end

opts = optimset('Display','iter','GradObj','on','TolX',10e-06);
A = [];
B = [];
Aeq = []; % ignore constant for linear constrant
beq = [];
lb = []; % ignore constant for positivity constraint
ub = [];
nonlincon = @(x) fitStim_nonlincon(x);
prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);

%-------------------------------------------
function [err grad] = fitStim_errFun(X_c,fit1,fit2)

% eval the output for model 1
[Y_st_1 dY_st_1 ddY_st_1] = evalNonLin((fit1.I_sc.*fit1.A_sc)*X_c,fit1.f);
[Z_t_1 dZ_t_1 ddZ_t_1] = evalNonLin((fit1.J_s.*fit1.B_s)*Y_st_1,fit1.g);

% eval the output for model 2
[Y_st_2 dY_st_2 ddY_st_2] = evalNonLin((fit2.I_sc.*fit2.A_sc)*X_c,fit2.f);
[Z_t_2 dZ_t_2 ddZ_t_2] = evalNonLin((fit2.J_s.*fit2.B_s)*Y_st_2,fit2.g);

err = ((Z_t_1 - Z_t_2)^2);

for ic=1:fit1.C
    grad1(ic) = dZ_t_1*(fit1.B_s*(dY_st_1.*fit1.A_sc(:,ic).*fit1.I_sc(:,ic)));
    grad2(ic) = dZ_t_2*(fit2.B_s*(dY_st_2.*fit2.A_sc(:,ic).*fit2.I_sc(:,ic)));
end

grad = 2*(Z_t_1-Z_t_2)*(grad1-grad2);

%-------------------------------------------
function [C Ceq] = fitStim_nonlincon(x)

C = [];
Ceq = sum(x.^2)-1;