function fit = fitF(data,fit,initFlag)

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

displayMode = 'off';

% evaluate input to subunits
Y_st_in = (fit.I_sc.*fit.A_sc)*data.X_ct;

% get the range of values for setting the spline knots
mnVal = prctile(Y_st_in(:),1);
mxVal = prctile(Y_st_in(:),99);

% store the old nonlinearity
fOld = fit.f;

% create the spline parameter matrix
fit.f.knots = linspace(mnVal,mxVal,fit.f.nknots);
fit.f.Mspline = splineParamMatrix(fit.f.knots,...
    fit.f.smoothness,fit.f.extrap);

% get initial values for the spline parameters
initPrs = initNonLinSpline(fit.f,fOld);

% do the optimization
if ~exist('initFlag','var') || initFlag == 0
    options = optimset('maxiter',250,'maxfunevals',1e6,...
        'Display',displayMode,'Largescale','off','TolX',10e-04);
    estParams = fminunc(@(prs) fitF_errFun(prs,Y_st_in,data.R_t,fit),initPrs,options);
    fit.f.w = estParams;
else
    fit.f.w = initPrs;
end


%-------------------------------------------
function [err] = fitF_errFun(prs,Y_st_in,R_t,fit)

fit.f.w = prs;

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
    case 'mse'
        err = sum((Z_t-R_t).^2);

end
        