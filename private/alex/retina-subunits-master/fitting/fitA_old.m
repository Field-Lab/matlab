function fit = fitA(data,fit,testMode)

%-------------------------------------------
%
% fit = fitA(data,fit,testMode)
%
% fit parameters of a LN cascade model
% specifically, fit the weights from stimuli to subunits,
% and any parameters of the subunit nonlinearity
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
% this function fits the weights A_sc and the parameters
% of the subunit nonlinearity f assuming all other parameters
% are fixed
%
% if 'testMode' == 1, the function just tests
% the derivitives and hessians and doesn't do
% the optimization
%
% freeman, 3-19-2012
%-------------------------------------------

if nargin < 3
    testMode = 0;
end

if nargin < 2
    error('(fitA) need to provide data and fit structs');
end

% initialize the parameters
initPrs = [vector(fit.A_sc(find(fit.I_sc))); fit.f.p(1)];
%initPrs = [fit.A_sc(find(fit.I_sc))];

% create handle to objective function
fun = @(prs) fitA_errFun(prs(:),data.X_ct,data.R_t,fit);
x0 = initPrs(:);

if testMode == 1
    DerivCheck_Elts(fun,x0);
    HessCheck(fun,x0);
    return
end

if strcmp(fit.constraints.A_sc,'L1')
    hessfun = @(prs,lambda) getOutputs(@() fitA_errFun(prs(:),data.X_ct,data.R_t,fit),3);
    opts = optimset('Display','iter','GradObj','on','Hessian','user-supplied','HessFcn',hessfun,'Algorithm','interior-point','TolX',10e-06);
    A = [];
    B = [];
    Aeq = buildConstraintL1(fit.I_sc);
    Aeq(end+1,end+1) = 0; % ignore constant for linear constrant
    beq = ones(1,fit.S);
    beq(end+1) = 0;
    lb = zeros(1,length(x0));
    lb(end) = -inf; % ignore constant for positivity constraint
    ub = [];
    nonlincon = [];
    prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);
end
if strcmp(fit.constraints.A_sc,'pos')
    hessfun = @(prs,lambda) getOutputs(@() fitA_errFun(prs(:),data.X_ct,data.R_t,fit),3);
    opts = optimset('Display','iter','GradObj','on','Hessian','on','Algorithm','trust-region-reflective','TolX',10e-06);
    A = [];
    B = [];
    Aeq = [];
    beq = [];
    lb = zeros(1,length(x0));
    lb(end) = -inf; % ignore constant for positivity constraint
    ub = [];
    nonlincon = [];
    prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);
end
if strcmp(fit.constraints.A_sc,'none')
    opts = optimset('Display','iter','GradObj','on','Hessian','on');
    prsEst = fminunc(fun,x0,opts);
end

% collect the estimated parameters
inds = find(fit.I_sc);
nPrs = sum(fit.I_sc(:));
fit.A_sc(inds) = prsEst(1:nPrs);
fit.f.p(1) = prsEst(nPrs+1:end);


%-------------------------------------------
function [err grad hess] = fitA_errFun(prs,X_ct,R_t,fit)

% extract and store the weights where indices are set to 1
prinds = fit.I_sc(:) == 1;
nPrs = sum(fit.I_sc(:));
fit.A_sc(prinds) = prs(1:nPrs);

% extract the nonlinearity parameters
fit.f.p(1:(length(prs)-nPrs)) = prs(nPrs+1:end);

% get the current subunit outputs
[Y_st dY_st ddY_st] = evalNonLin((fit.I_sc.*fit.A_sc)*X_ct,fit.f);

% get the current RGC output
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
        
        if nargout > 1
            % get gradient
            W_s = fit.J_s.*fit.B_s;
            dY_st_W_s = bsxfun(@times,dY_st,W_s');
            grad_A = -(bsxfun(@times,dY_st_W_s,R_t-dZ_t)*X_ct');
            grad_A = grad_A(prinds);
            grad_mu = -sum((Z_t-R_t).*(W_s*dY_st));
            grad = [vector(grad_A); grad_mu];
        end
        
        if nargout > 2
            % get hessian
            i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
            j = 1:size(X_ct,2); j = j(ones(1,1),:);
            mult1 = X_ct(i,j);
            mult2 = repmat(dY_st_W_s,[fit.C 1]);
            all_vec_1 = mult1(prinds,:) .* mult2(prinds,:);
            H_A = -all_vec_1*bsxfun(@times,all_vec_1,-Z_t)';
            H_mu = sum((Z_t).*((W_s*dY_st).^2));
            H_mu_A = -sum(bsxfun(@times,all_vec_1,(Z_t).*(W_s*dY_st)),2);
            hess = [H_A H_mu_A; H_mu_A', H_mu];
            %hess = H_A;           
        end
        
    otherwise
        if nargout > 1
            % get gradient
            Q_t = dZ_t./Z_t;
            W_s = fit.J_s.*fit.B_s;
            dY_st_W_s = bsxfun(@times,dY_st,W_s');
            grad_A = -(bsxfun(@times,dY_st_W_s,(R_t.*Q_t)-dZ_t)*X_ct');
            grad_A = grad_A(prinds);
            grad_mu = -sum((dZ_t-R_t.*(dZ_t./Z_t)).*(W_s*dY_st));
            grad = [vector(grad_A); grad_mu];
        end
        
        if nargout > 2
            % get hessian
            Q_t = ((Z_t.*ddZ_t)-(dZ_t.*dZ_t))./(Z_t.*Z_t);
            W_s = fit.J_s.*fit.B_s;
            dY_st_W_s = bsxfun(@times,dY_st,W_s');
            i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
            j = 1:size(X_ct,2); j = j(ones(1,1),:);
            mult = X_ct(i,j);
            all_vec_1 = repmat(dY_st_W_s,[fit.C 1]) .* mult;
            H_B_1 = all_vec_1*bsxfun(@times,all_vec_1,(R_t.*Q_t)-ddZ_t)';
            tmpinds = find(fit.I_sc);
            hess = -H_B_1(tmpinds,tmpinds);
            % this part goes away if f is rectified linear
            %Q_t = (dZ_t./Z_t)
            %ddY_st_W_s = bsxfun(@times,ddY_st,W_s');
            %i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
            %j = 1:size(X_ct,2); j = j(ones(1,1),:);
            %mult = X_ct(i,j);
            %mul2 = repmat(ddY_st_W_s,[fit.C 1])
            %H_B_2 = mult*bsxfun(@times,mult.*mult2,((R_t.*Q_t)-dZ_t))';
            
        end
end