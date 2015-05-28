function [muc Sigmac sample]=SampleConditionalNormal(mu,Sigma,x,ind2)

ind1=setdiff([1:length(mu)],ind2);
mu1=mu(ind1);
mu2=mu(ind2);
Sigma11=Sigma(ind1,ind1);
Sigma22=Sigma(ind2,ind2);
Sigma12=Sigma(ind1,ind2);

muc=mu2+Sigma12'*Sigma11^(-1)*(x-mu1);

if(nargout>1)
Sigmac=Sigma22-Sigma12'*Sigma11^-1*Sigma12;
Sigmac=(Sigmac+Sigmac')/2;
sample=mvnrnd(muc',full(Sigmac))';
end
