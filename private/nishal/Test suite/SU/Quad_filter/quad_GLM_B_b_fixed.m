
function [f,g]=quad_GLM_B_b_fixed(x)
global maskedMov binnedResponses masksz b rankB

c=x(1);
%b=x(2:2+masksz-1);
b_sc = x(2);
B=reshape(x(3:end),[masksz,rankB]);

data_short=maskedMov';
dt=1/120;
resp=binnedResponses;
n=length(resp);

lam = data_short*b*b_sc + c + sum((data_short*B).^2,2);
lam=lam';
lam_exp=exp(0.15*lam);

%% function value
f=+(sum(lam_exp*dt) - sum(0.15*lam*resp));

%% gradients
g(1) = sum(lam_exp*dt) - sum(resp);

gb = sum((data_short*b).*repmat(lam_exp'*dt - resp,[1,1]),1); 
g(2) = gb;

gB = ((repmat ((-resp + lam_exp'*dt),[1,masksz]).*data_short)'*data_short)*2*B;
g(3:3+ masksz*rankB-1) = gB(:);

end

