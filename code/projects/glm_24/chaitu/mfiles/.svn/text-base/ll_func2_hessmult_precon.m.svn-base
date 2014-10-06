% It returns function, gradient (and Hessian?) value of the log likelihood
% of the data. This function can be used by fmincon

function [f g Hinfo] = ll_func2_hessmult_precon(p,basepars,stimpars,trainpars,Mprecon)
1;
% Calculate the log likelihood of the data
[logprob cifs] = train_ll3(Mprecon*p,basepars,stimpars,trainpars);
f = full(-logprob); % remember to negate - we want to maximize log likelihood
% Calculate gradient
if (nargout < 3)
    [g] = train_ll_grad4(Mprecon*p,basepars,stimpars,trainpars,cifs);
else
    [g Hinfo] = train_ll_grad_hessmult(Mprecon*p,basepars,stimpars,trainpars,cifs);
end
g = -Mprecon'*g;
1;
fprintf('f=%f    gnorm=%f\n',f,norm(g));
