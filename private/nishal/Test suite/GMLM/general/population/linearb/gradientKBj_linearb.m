function [gradK,gradB] = gradientKBj_linearb(X,Y,Bs,Ks,T,Nc,alphak,alphab,dt)


sta = zeros(size(X,1),1);
for icell=1:Nc
 sta = sta + X*(alphak(:,icell)'.*Y(icell,:))'/T;
end

const = sum(Bs)*dt/T;



gradK = const*X*exp(Ks'*X)' - sta;
gradB = (dt/T)*sum(exp(Ks'*X))- sum(alphab'.*Y,2)'/T;


end