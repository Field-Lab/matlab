% It returns function, gradient (and Hessian) value of the log likelihood
% of the data with respect to the base firing rate, 
% two independent magnitude coefficients for two filter kernels (taken from rank-2 GLM or STA),
% and postspike/coupling filters coeffs.
% This function can be used by fmincon

function [f g H] = ll_func2_fixrk2(p,basepars,stimpars,trainpars)
% Calculate the log likelihood of the data

[logprob lcifs cifs] = train_ll3_fixrk2(p,basepars,stimpars,trainpars);
f = full(-logprob); % remember to negate - we want to maximize log likelihood
% Calculate gradient
g = zeros(size(p));
kernel1 = full((1.*trainpars.D) - trainpars.dt.*cifs);

% Calculate gradient
if (nargout < 3)
    [g] = train_ll_grad4_fixrk2(p,basepars,stimpars,trainpars,lcifs,cifs);
else
    [g H] = train_ll_grad4_fixrk2(p,basepars,stimpars,trainpars,lcifs,cifs);
    H = -H;
    %fprintf('Condition number of total hessian is %f\n',cond(H));
end
g = -g;

fprintf('f=%f    gnorm=%f\n',f,norm(g));
1;
