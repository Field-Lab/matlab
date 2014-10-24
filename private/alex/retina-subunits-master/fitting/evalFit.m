function [out] = evalFit(data,fit,inds)

%-------------------------------------------
%
% evalFit(data,fit)
%
% get various statistics of the model fit
%
% inputs:
% data -- stucture containing the data
% fit -- structure containing fit parameters
%
% outputs:
% out -- structure with the following fields
%        'negloglik' -- negative log likelihood
%        'r2' -- r squared
%        'Z_t' -- predicted response
%
% freeman, 1-1-2012
%-------------------------------------------

if ~exist('inds','var') || isempty(inds)
    inds = 1:length(data.R_t);
end

switch fit.model
    case 'subunit'
        % get the subunit responses
        [Y_st] = evalNonLin((fit.I_sc.*fit.A_sc)*data.X_ct,fit.f);
        % get the current RGC output
        [Z_t] = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);
        
    case 'stc'
        % get the STA filter response
        K_st = fit.k'*data.X_ct;
        % get the STC filter responses
        V_st = (fit.V(:,1:fit.D-1)'*data.X_ct).^2;
        % combine
        input_st = [K_st;V_st];
        % get weighted sum of STA and STC responses
        Z_t = evalNonLin(fit.d_s*input_st,fit.g);
    case 'ln'
        Z_t = evalNonLin(fit.d_s*data.X_ct,fit.g);
        
end

% check for small values, set to 0 for lik
etol = 10^-6;
iiz = Z_t <= etol;
Z_t(iiz) = etol;

% compute liklihood
loglik = getLogLikSpk(Z_t(inds),data.R_t(inds));
negloglik = -loglik;

% compute r2
sse = sum((Z_t(inds)-data.R_t(inds)).^2);
sst = sum((data.R_t(inds)-mean(data.R_t(inds))).^2);
r2 = 1 - sse/sst;

out.negloglik = negloglik;
out.r2 = r2;
out.Z_t = Z_t;

