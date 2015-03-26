function [f g] = test_sli_hess(x,pars)

[realf g Hinfo] = sli_wrapper(x,pars);
f = g(pars.hesstest_j);
ej = zeros(size(x)); ej(pars.hesstest_j) = 1;
g = sli_newton_hessmult2(Hinfo,ej,pars.Kfunc,pars.E,pars);
