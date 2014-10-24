function fit = fitLN(data,fit,testMode)

if nargin < 3
    testMode = 0;
end

if nargin < 2
    error('(fitLN) need to provide data and fit structs');
end

% collect filter outputs
input_st = data.X_ct;
nPrs = size(input_st,1);

% initialize the parameters
initPrs = [fit.d_s fit.g.p(1)];

fun = @(prs) fitLN_errFun(prs(:),input_st,data.R_t,fit);
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
fit.d_s(1:nPrs) = prsEst(1:nPrs);

% collect the output nonlinearity parameters
fit.g.p(1) = prsEst(nPrs+1:end);

%-------------------------------------------
function [err grad] = fitLN_errFun(prs,input_st,R_t,fit)

% get number of inputs
nPrs = size(input_st,1);

% grab current value of B_s
fit.d_s(1:nPrs) = prs(1:nPrs);

% grab curent nonlinearity parameters
fit.g.p(1:(length(prs)-nPrs)) = prs(nPrs+1:end);


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


