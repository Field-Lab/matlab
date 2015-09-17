function [gradK,gradB] = gradientKBj(X,Y,Bs,Ks,T,Nc,alpha,dt)


sta = zeros(size(X,1),1);
for icell=1:Nc
    a_use = reshape(alpha(1,:,icell),[T,1]);
 sta = sta + X*(a_use'.*Y(icell,:))'/T;
end

const = sum(exp(Bs))*dt/T;



gradK = const*X*exp(Ks'*X)' - sta;
a = alpha(1,:,:);
a =permute(a,[2,3,1]);
gradB = exp(Bs)*(dt/T)*sum(exp(Ks'*X))- sum(a'.*Y,2)'/T;


end