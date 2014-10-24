function fit = fitSTC(data,fit,testMode)

%-------------------------------------------
%
% fit = fitSTC(data,fit,testMode)
%
% fit a weighted sum of STA and STC components
% to preidct the response of a cell
%
% data should contain the fields:
% data.X_ct -- stimulus, C x T
% data.R_t -- stimulus, 1 x T
%
% fit should contain the fields:
% fit.k -- sta
% fit.V -- matrix of stc filters
% fit.e -- eigenvalues
%
% this function fits a set of weights on the STA
% STC filters, along with a constant, to maximize
% the likelihood of the ganglion cell response
%
% if 'testMode' == 1, the function just tests
% the derivitives and hessians and doesn't do
% the optimization
%
% freeman, 4-1-2012
%-------------------------------------------


if nargin < 3
    testMode = 0;
end

if nargin < 2
    error('(fitSTC) need to provide data and fit structs');
end

% collect filter outputs
K_st = fit.k'*data.X_ct;
if fit.D > 1
    V_st = (fit.V(:,1:fit.D-1)'*data.X_ct).^2;
else
    V_st = [];
end
input_st = [K_st;V_st];

% initialize the parameters
initPrs = [fit.d_s fit.g.p(1)];

% create handle to objective function
fun = @(prs) fitSTC_errFun(prs(:),input_st,data.R_t,fit);
x0 = initPrs(:);

% check gradients
if testMode == 1
    DerivCheck_Elts(fun,x0);
    %HessCheck(fun,x0);
    return
end

opts = optimset('Display','off','GradObj','on');
prsEst = fminunc(fun,x0,opts);

% collect the estimated B_s parameters
fit.d_s(1:fit.D) = prsEst(1:fit.D);

% collect the output nonlinearity parameters
fit.g.p(1) = prsEst(fit.D+1:end);

%-------------------------------------------
function [err grad] = fitSTC_errFun(prs,input_st,R_t,fit)


% grab current value of B_s
fit.d_s(1:fit.D) = prs(1:fit.D);

% grab curent nonlinearity parameters
fit.g.p(1:(length(prs)-fit.D)) = prs(fit.D+1:end);

% get the predicted response given these parameters
[Z_t dZ_t ddZ_t] = evalNonLin((fit.d_s*input_st),fit.g);

% check for small values, set to 0 for lik and gradient
etol = 10^-6;
iiz = Z_t <= etol;
Z_t(iiz) = etol;
dZ_t(iiz) = 0;
ddZ_t(iiz) = 0;

% compute liklihood
loglik = getLogLikSpk(Z_t,R_t);
err = -loglik;

switch fit.g.type
    case 'exp' % simpler form when g is exponential
        
        % compute gradient
        if nargout > 1
            grad_gamma = -(sum(Z_t-R_t)); % constant grad
            grad_d = -((R_t-Z_t)*input_st'); % weight grad
            grad = [grad_d grad_gamma]; % combined grad
        end
    otherwise
        % compute gradient
        if nargout > 1
            Q_t = dZ_t./Z_t; % common terms
            grad_gamma = -(sum(dZ_t-(R_t.*Q_t))); % constant grad
            grad_d = -(((R_t.*Q_t)-dZ_t)*input_st'); % weight grad
            grad = [grad_d grad_gamma]; % combined grad
        end
end


