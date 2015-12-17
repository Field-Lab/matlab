function [pred,lam]=response_lnp(data,stimulus,nTrials)

 mov=stimulus; % 3D, between -0.5 to +0.5.
maskedMovdd= filterMov(mov,data.mask,squeeze(data.ttf));
inputt = (data.usta'*maskedMovdd - data.ainput) / data.binput;

frnz = data.fr >0; 
f=fit(data.in(frnz),log(data.fr(frnz)),'poly2');
g = @(x) exp(f(x));
%  subplot(1,2,2);
%  plot(data.in,data.fr,data.in,g(data.in)); pause(0.1);
 
lam = min(g(inputt),max(data.fr));
lam = lam*data.dt;


pred=zeros(nTrials,length(lam));
for itrial=1:nTrials
pred(itrial,:)=poissrnd(lam);    
end

 
end