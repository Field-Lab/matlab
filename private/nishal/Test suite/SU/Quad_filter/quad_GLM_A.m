
function [f,g]=quad_GLM_A(x)
global maskedMov binnedResponses masksz gamma wt

c=x(1);
b=x(2:2+masksz-1);
A=reshape(x(2+masksz:end),[masksz,masksz]);
wt_sq = reshape(wt,[masksz,masksz]);

data_short=maskedMov';
dt=1/120;
resp=binnedResponses;
n=length(resp);


lam = data_short*b + c + sum((A*data_short').*data_short',1)';
lam=lam';


lam_exp=10*exp(0.15*lam);

%% function value

f=+(sum(lam_exp*dt) - sum(0.15*lam*resp)) + sum(gamma*abs(wt.*A(:)));

%% gradients

g(1) = sum(lam_exp*dt) - sum(resp);

gb = sum(data_short.*repmat(lam_exp'*dt - resp,[1,masksz]),1); 
g(2:2+masksz-1) = gb;

gA= ((repmat ((-resp + lam_exp'*dt),[1,masksz]).*data_short)'*data_short) + gamma*wt_sq.*sign(A);

g(2+masksz:2+masksz + masksz^2-1) = gA(:);

end
