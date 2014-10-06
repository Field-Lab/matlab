function y=sinus_fit(beta,t);
a=size(beta);
if a(1)>1
    beta=beta';
end
y=beta(1,1)*sin(2*pi*beta(1,2)*t+beta(1,3))+beta(1,4);