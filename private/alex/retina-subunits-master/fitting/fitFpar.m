function fit = fitFpar(data,fit,initFlag)

%-------------------------------------------
%
% fit = fitF(data,fit,testMode)
%
% fit parameters of a LN cascade model
% specifically, fit the subunit nonlinearity
%
% data should contain the fields:
% data.X_ct -- stimulus, C x T
% data.R_t -- stimulus, 1 x T
%
% fit should contain the fields:
% fit.S -- number of subunits
% fit.I_sc -- binary assignments from stimulus to subunits, S x C
% fit.A_sc -- weights from stimulus to subunits, S x C
% fit.J_s -- binary assignments from subunits to cell, 1 x S
% fit.B_s -- weights from subunits to cell, 1 x S
% fit.f -- subunit nonlinearity (see evalNonLin.m)
% fit.g -- output nonlinearity (see evalNonLin.m)
%
% this function fits the fit.f as a cubic spline
%
% freeman, 3-14-2012
%-------------------------------------------

displayMode = 'iter';

% evaluate input to subunits
Y_st_in = (fit.I_sc.*fit.A_sc)*data.X_ct;
opts = optimset('Display',displayMode,'MaxIter',100,'maxfunevals',1e6,'GradObj','off');
fun = @(prs) fitFpar_errFun(prs,Y_st_in,data.R_t,fit);
x0 = fit.f.p;
prsEst = fminunc(fun,x0,opts);
fit.f.p = prsEst;

%-------------------------------------------
function [err] = fitFpar_errFun(prs,Y_st_in,R_t,fit)

fit.f.p = prs;

% get the current subunit outputs
Y_st = evalNonLin(Y_st_in,fit.f);

% get the current RGC output
Z_t = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);

switch fit.error
    case 'loglik'
        % check for small values, set to 0 for gradient
        etol = 10^-6;
        iiz = Z_t <= etol;
        f(iiz) = etol;
        df(iiz) = 0;
        dff(iiz) = 0;
        
        % compute liklihood
        loglik = getLogLikSpk(Z_t,R_t);
        err = -loglik;

end
