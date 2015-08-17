function [LL,grad,Hess]=fit_divgainASM2_grad_gain_onlynr_cvx(y,X,Y,params,su,f,fd,fd2,interval,staf,dnr)


isu=su;
params2.k.value = y(params{isu}.k.idx2);


%% compute likelihood
T=length(Y);
innerprod = params2.k.value'*X;
total_inp2 =  f((innerprod)) ./dnr;    
LL = sum(total_inp2)*interval/(120*T) - params2.k.value'*staf;

%% compute gradient

    fdd = fd(innerprod)./dnr;
    params2.k.grad  = X*fdd'*interval/(120*T) - staf;

 grad = gpuArray(0*y);
 grad(params{isu}.k.idx2) = params2.k.grad;

 hessnr = fd2(innerprod)./dnr;
 Hess = (interval/(120*T))*(repmat(hessnr,[size(X,1),1]).*X)*X';

LL = gather(LL);
grad =gather(grad);
Hess = gather(Hess);
end