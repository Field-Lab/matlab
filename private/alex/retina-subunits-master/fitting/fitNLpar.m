function nl = fitNLpar(x,y,nl)

displayMode = 'off';

x0 = rand(5,1);
%x0 = [-1.8 -2 -0.3 -0.3];

opts = optimset('Display',displayMode,'MaxFunEvals',10E6,'GradObj','off');
fun = @(prs) fitNLpar_errFun(prs,x(:),y(:),nl);

%keyboard
%DerivCheck_Elts(fun,x0);
[prsEst,err] = fminunc(fun,x0,opts);

nl.p = prsEst;

function [sse grad] = fitNLpar_errFun(prs,x,y,nl)

nl.p = prs;
yhat = evalNonLin(x,nl);
sse = sum((yhat-y).^2);

if nargout > 1
    a = prs(1);
    b = prs(2);
    d = prs(3);
    e = prs(4);
    %e = prs(5);
    diff = y-yhat;
    eterm = exp(b*x+a);
    grad_a_tmp = exp(d)*(log(1+eterm)).^(exp(d)-1) .* (eterm)./(1+eterm);
    grad_a = -sum(2*(diff).*grad_a_tmp);
    grad_b_tmp = exp(d)*(log(1+eterm)).^(exp(d)-1) .* (x.*eterm)./(1+eterm);
    grad_b = -sum(2*(diff).*grad_b_tmp);
    %grad_c_tmp = (log(1+exp(b*x+a))).^d;
    %grad_c = -sum(2*(diff).*grad_c_tmp);
    grad_d_tmp = (log(1+eterm)).^exp(d) .* log(log(1+eterm)) * exp(d);
    grad_d = -sum(2*(diff).*grad_d_tmp);
    grad_e_tmp = 1;
    grad_e = -sum(2*(diff).*grad_e_tmp);
    grad = [grad_a grad_b grad_d grad_e];
end
    
