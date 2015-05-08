
function [f,g]=quad_GLM_A_b_fixed(x)
global maskedMov binnedResponses masksz gamma wt b_ratios

c=x(1);
b_sc=x(2);

b=b_ratios*b_sc;

A=reshape(x(3:end),[masksz,masksz]);
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

gb = sum(sum(data_short.*repmat(lam_exp'*dt - resp,[1,masksz]),1)); 
g(2) = gb;

gA= ((repmat ((-resp + lam_exp'*dt),[1,masksz]).*data_short)'*data_short) + gamma*wt_sq.*sign(A);

g(3:3 + masksz^2-1) = gA(:);

end
