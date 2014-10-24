function fit = fitG(data,fit,initFlag)

%-------------------------------------------
%
% fit = fitG(data,fit,initFlag)
%
% fit parameters of a LN cascade model
% specifically, fit the output nonlinearity
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
% this function fits the fit.g as a cubic spline
%
% if 'testMode' == 1, the function just tests
% the derivitives and hessians and doesn't do
% the optimization
%
% freeman, 3-14-2012
%-------------------------------------------

switch fit.model
     case 'subunit'
        % evaluate input to RGC under subunit model
        Y_st = evalNonLin((fit.I_sc.*fit.A_sc)*data.X_ct,fit.f);
        Z_st_in = (fit.J_s.*fit.B_s)*Y_st;
        
     case 'stc'
        % evaluate input to RGC under stc model
        K_st = fit.k'*data.X_ct;
        if fit.D > 1
            V_st = (fit.V(:,1:fit.D-1)'*data.X_ct).^2;
        else
            V_st = [];
        end
        Z_st_in = (fit.d_s*[K_st;V_st]);
        
    case 'ln'
        % evaluate input to RGC under stc model
        Z_st_in = (fit.d_s*data.X_ct);
end

% get the range of values for setting the spline knots
mnVal = prctile(Z_st_in(:),1);
mxVal = prctile(Z_st_in(:),99);

% store the old nonlinearity
fOld = fit.g;

% create the spline parameter matrix
fit.g.knots = linspace(mnVal,mxVal,fit.g.nknots);
fit.g.Mspline = splineParamMatrix(fit.g.knots,...
    fit.g.smoothness,fit.g.extrap);

% get initial values for the spline parameters
initPrs = initNonLinSpline(fit.g,fOld);

% do the optimization
if ~exist('initFlag','var') || initFlag == 0
    options = optimset('maxiter',250,'maxfunevals',1e6,...
        'Display','off','Largescale','off');
    estParams = fminunc(@(prs) fitF_errFun(prs,Z_st_in,data.R_t,fit),initPrs,options);
    fit.g.w = estParams;
else
    fit.g.w = initPrs;
end


%-------------------------------------------
function [err] = fitF_errFun(prs,Z_st_in,R_t,fit)

fit.g.w = prs;

% get the current RGC output
Z_t = evalNonLin(Z_st_in,fit.g);

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
        