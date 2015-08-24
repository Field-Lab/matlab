function [LL,grad,Hess] = optmizeK(k,X,const,sta)
LL = const*sum(exp(k'*X)) - k'*sta;

grad = const*X*exp(k'*X)' - sta; 

Hess = const*((repmat(exp(k'*X),[size(X,1),1]).*X)*X');

LL = gather(LL);
grad = gather(grad);
Hess = gather(Hess);
end