function fit = fitB(data,fit,testMode)

%-------------------------------------------
%
% fit = fitB(data,fit,testMode)
%
% fit parameters of a LN cascade model
% specifically, fit the weights from subunits to cell,
% and any parameters of the output nonlinearity
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
% this function fits the weights B_s and the parameters
% of the output nonlinearity g assuming all other parameters
% are fixed
%
% if 'testMode' == 1, the function just tests
% the derivitives and hessians and doesn't do
% the optimization
%
% freeman, 3-14-2012
%-------------------------------------------

displayMode = 'off';

if nargin < 3
    testMode = 0;
end

if nargin < 2
    error('(fitB) need to provide data and fit structs');
end

% evaluate the subunit responses
Y_st = evalNonLin((fit.I_sc.*fit.A_sc)*data.X_ct,fit.f);

if strcmp(fit.g.type,'spline')
    usingSplines = 1;
else
    usingSplines = 0;
end

% initialize the parameters
if usingSplines
    initPrs = [fit.B_s];
else
    initPrs = [fit.B_s, fit.g.p(1)];
end

% create handle to objective function
fun = @(prs) fitB_errFun(prs(:),Y_st,data.R_t,fit);
x0 = initPrs(:);

if testMode == 1
    DerivCheck_Elts(fun,x0);
    HessCheck(fun,x0);
    return
end

% set up the optimization, depends on the constraints
if strcmp(fit.constraints.B_s,'L1')
    hessfun = @(prs,lambda) getOutputs(@() fitB_errFun(prs(:),Y_st,data.R_t,fit),3);
    opts = optimset('Display',displayMode,'GradObj','on','Hessian','user-supplied','HessFcn',hessfun,'Algorithm','interior-point');
    A = [];
    B = [];
    Aeq = [];
    beq = [];
    %Aeq = ones(1,length(x0));
    %Aeq(end,end) = 0; % ignore constant for linear constraint
    %beq = 1;
    lb = zeros(1,length(x0));
    lb(end) = -inf; % ignore constant for positivity constraint
    ub = [];
    nonlincon = [];
    prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);
end
if strcmp(fit.constraints.B_s,'none')
    opts = optimset('Display',displayMode,'GradObj','on','Hessian','on');
    prsEst = fminunc(fun,x0,opts);
end

% collect the estimated B_s parameters
fit.B_s(1:fit.S) = prsEst(1:fit.S);

if ~usingSplines
    % collect the output nonlinearity parameters
    fit.g.p(1) = prsEst(fit.S+1:end);
end

tmpOut = (fit.J_s.*fit.B_s)*evalNonLin((fit.I_sc.*fit.A_sc)*data.X_ct,fit.f);
fit.g.rng = [min(tmpOut(:)) max(tmpOut(:))];

%-------------------------------------------
function [err grad hess] = fitB_errFun(prs,Y_st,R_t,fit)

if strcmp(fit.g.type,'spline')
    usingSplines = 1;
else
    usingSplines = 0;
end

% grab current value of B_s
fit.B_s(1:fit.S) = prs(1:fit.S);

% grab curent nonlinearity parameters
if ~usingSplines
    fit.g.p(1:(length(prs)-fit.S)) = prs(fit.S+1:end);
end

% get the predicted response given these parameters
[Z_t dZ_t ddZ_t] = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);

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
            grad_B = -((R_t-Z_t)*Y_st'); % weight grad
            if usingSplines
                grad = [grad_B];
            else
                grad_gamma = -(sum(Z_t-R_t)); % constant grad
                grad = [grad_B grad_gamma]; % combined grad
            end
        end
        
        % compute hessian
        if nargout > 2  
            % common terms
            Y_st_J_s = bsxfun(@times,Y_st,fit.J_s');
            % weight hess
            H_B = -(Y_st_J_s*bsxfun(@times,Y_st_J_s,-Z_t)');
            % combined hess
            if usingSplines
                hess = [H_B];
            else
                % constant hess
                H_gamma = -(sum(-Z_t));
                H_gamma_B = -((Z_t)*Y_st_J_s');
                hess = [H_B H_gamma_B'; H_gamma_B H_gamma];
            end
        end 
        
    otherwise % this is for a generic g
        
        % compute gradient
        if nargout > 1
            Q_t = dZ_t./Z_t; % common terms
            grad_B = -(((R_t.*Q_t)-dZ_t)*Y_st'); % weight grad
            if usingSplines
                grad = [grad_B];
            else
                grad_gamma = -(sum(dZ_t-(R_t.*Q_t))); % constant grad
                grad = [grad_B grad_gamma]; % combined grad
            end
        end
        
        % compute hessian
        if nargout > 2  
            % common terms
            Y_st_J_s = bsxfun(@times,Y_st,fit.J_s');
            Q_t = (Z_t.*ddZ_t-dZ_t.*dZ_t)./(Z_t.*Z_t);
            % weight hess
            H_B = -(Y_st_J_s*bsxfun(@times,Y_st_J_s,(R_t.*Q_t)-ddZ_t)');
            if usingSplines
                hess = [H_B];
            else
                % constant hess
                H_gamma = -(sum((R_t.*Q_t) - ddZ_t));
                H_gamma_B = -((ddZ_t-(R_t.*Q_t))*Y_st_J_s');
                % combined hess
                hess = [H_B H_gamma_B'; H_gamma_B H_gamma];
            end
        end 
end