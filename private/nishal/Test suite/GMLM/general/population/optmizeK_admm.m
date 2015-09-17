
function [LL,grad,Hess] = optmizeK_admm(x,X,sta,b_coeff,rho,K_zk,u_Kk,Kidx,Bidx,dt,scale_cell)

k = x(Kidx);
b = x(Bidx);

T =size(X,2);

const = sum(exp(b).*(T./scale_cell))*dt/T;

% Get LL
LL = const*sum(exp(k'*X)) - k'*sta -b_coeff'*b + (rho/2)*norm(k - K_zk + u_Kk)^2;

% Get grad
grad = gpuArray(zeros(length(x),1));
grad(Kidx) = const*X*exp(k'*X)' - sta +rho*(k - K_zk + u_Kk); 
grad(Bidx) = (exp(b).*(T./scale_cell))*(dt/T)*sum(exp(k'*X)) - (b_coeff);

% Get Hessian
Hess = gpuArray(zeros(length(x),length(x)));
Hess(Kidx,Kidx) = const*((repmat(exp(k'*X),[size(X,1),1]).*X)*X') + (rho)*eye(length(k),length(k));
Hess(Bidx,Bidx) = diag((exp(b).*(T./scale_cell))*(dt/T)*sum(exp(k'*X)));
Hess(Kidx,Bidx) = (dt/T)*X*exp(k'*X)'*(exp(b).*(T./scale_cell))';
Hess(Bidx,Kidx) = Hess(Kidx,Bidx)';

LL = gather(LL);
grad = gather(grad);
Hess = gather(Hess);
end