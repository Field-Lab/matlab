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

displayMode = 'off';

if nargin < 3
    testMode = 0;
end

if nargin < 2
    error('(fitA) need to provide data and fit structs');
end

% determine whether we are using splines
% (and don't need to estimate the constant)
if strcmp(fit.f.type,'spline') || strcmp(fit.f.type,'loglinear');
    usingSplines = 1;
else
    usingSplines = 0;
end

% get the indices of A_sc for which I_sc is 1
%prinds = fit.I_sc(:) == 1;
%nPrs = sum(fit.I_sc(:));
prinds = bsxfun(@and,fit.I_sc==1,sum(fit.I_sc,2)>1);
nPrs = sum(prinds(:));
inds_fixed = sum(prinds,2)==0;

% initialize the parameters
if usingSplines
    initPrs = [vector(fit.A_sc(prinds))];
else
    initPrs = [vector(fit.A_sc(prinds)); fit.f.p(1)];
end

% precompute subunit outputs that will stay fixed
[Y_st_fixed dY_st_fixed ddY_st_fixed] = evalNonLin((fit.I_sc(inds_fixed,:).*...
    fit.A_sc(inds_fixed,:))*data.X_ct,fit.f);

% create handle to objective function
fun = @(prs) fitA_errFun(prs(:),data.X_ct,data.R_t,fit,Y_st_fixed,dY_st_fixed,ddY_st_fixed);
x0 = initPrs(:);

if testMode == 1
    DerivCheck_Elts(fun,x0);
    HessCheck(fun,x0);
    %keyboard
    return
end

if strcmp(fit.constraints.A_sc,'L1')
    hessfun = @(prs,lambda) getOutputs(@() fitA_errFun(prs(:),data.X_ct,data.R_t,fit,Y_st_fixed,dY_st_fixed,ddY_st_fixed),3);
    if usingSplines
        %opts = optimset('Display','iter','GradObj','on','Hessian','user-supplied','HessFcn',hessfun,'Algorithm','interior-point','TolX',10e-04);
        opts = optimset('Display',displayMode,'GradObj','on','Algorithm','interior-point','TolX',10e-06);
    else
        opts = optimset('Display',displayMode,'GradObj','on','Hessian','user-supplied','HessFcn',hessfun,'Algorithm','interior-point','TolX',10e-06);
        %opts = optimset('Display','iter','GradObj','on','Algorithm','interior-point','TolX',10e-06);
    end
    A = [];
    B = [];
    Aeq = buildConstraintL1(fit.I_sc(sum(fit.I_sc,2)>1,:));
    beq = ones(1,sum(sum(fit.I_sc,2)>1));
    lb = zeros(1,length(x0));
    if ~usingSplines
        beq(end+1) = 0;
        Aeq(end+1,end+1) = 0; % ignore constant for linear constrant
        lb(end) = -inf; % ignore constant for positivity constraint
    end
    ub = [];
    nonlincon = [];
    prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);
end
if strcmp(fit.constraints.A_sc,'pos')
    hessfun = @(prs,lambda) getOutputs(@() fitA_errFun(prs(:),data.X_ct,data.R_t,fit),3);
    opts = optimset('Display',displayMode,'GradObj','on','Hessian','on','Algorithm','trust-region-reflective','TolX',10e-06);
    A = [];
    B = [];
    Aeq = [];
    beq = [];
    lb = zeros(1,length(x0));
    if ~strcmp(fit.f.type,'spline')
        lb(end) = -inf; % ignore constant for positivity constraint
    end
    ub = [];
    nonlincon = [];
    prsEst = fmincon(fun,x0,A,B,Aeq,beq,lb,ub,nonlincon,opts);
end
if strcmp(fit.constraints.A_sc,'none')
    if usingSplines
        opts = optimset('Display',displayMode,'GradObj','on','Hessian','off');
    else
        opts = optimset('Display',displayMode,'GradObj','on','Hessian','on');
    end
    prsEst = fminunc(fun,x0,opts);
end

% collect the estimated parameters
fit.A_sc(prinds) = prsEst(1:nPrs);

% grab parameters of nonlinearity 
% (unless we are using a spline)
if ~usingSplines
    fit.f.p(1) = prsEst(nPrs+1:end);
end

%-------------------------------------------
function [err grad hess] = fitA_errFun(prs,X_ct,R_t,fit,Y_st_fixed,dY_st_fixed,ddY_st_fixed)

% extract and store the weights where indices are set to 1
% and there are more than
%prinds = fit.I_sc(:) == 1;
%nPrs = sum(fit.I_sc(:));
prinds = bsxfun(@and,fit.I_sc==1,sum(fit.I_sc,2)>1);
nPrs = sum(prinds(:));
inds_fixed = sum(prinds,2)==0;

fit.A_sc(prinds) = prs(1:nPrs);

% determine whether we are using splines
% (and don't need to estimate the constant)
if (strcmp(fit.f.type,'spline') && ~isempty(fit.f.w)) || strcmp(fit.f.type,'loglinear')
    usingSplines = 1;
else
    usingSplines = 0;
end

% extract the nonlinearity parameters
if ~usingSplines
    fit.f.p(1:(length(prs)-nPrs)) = prs(nPrs+1:end);
end

if nargout == 1 % error only
    % get the current subunit outputs
    Y_st = zeros(fit.S,size(X_ct,2));
    [Y_st_sub] = evalNonLin((fit.I_sc(~inds_fixed,:).*fit.A_sc(~inds_fixed,:))*X_ct,fit.f);
    Y_st(inds_fixed,:) = Y_st_fixed;
    Y_st(~inds_fixed,:) = Y_st_sub;
    %[Y_st] = evalNonLin((fit.I_sc.*fit.A_sc)*X_ct,fit.f);
    % get the current RGC output
    [Z_t] = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);
end

if nargout == 2 % gradient only
    % get the current subunit outputs
    Y_st = zeros(fit.S,size(X_ct,2));
    dY_st_sub = Y_st;
    [Y_st_sub dY_st_sub] = evalNonLin((fit.I_sc(~inds_fixed,:).*fit.A_sc(~inds_fixed,:))*X_ct,fit.f);
    Y_st(inds_fixed,:) = Y_st_fixed;
    Y_st(~inds_fixed,:) = Y_st_sub;
    dY_st(inds_fixed,:) = dY_st_fixed;
    dY_st(~inds_fixed,:) = dY_st_sub;
    %[Y_st dY_st] = evalNonLin((fit.I_sc.*fit.A_sc)*X_ct,fit.f);
    % get the current RGC output
    [Z_t dZ_t] = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);
end

if nargout == 3 % gradient and hessian
    % get the current subunit outputs
    Y_st = zeros(fit.S,size(X_ct,2));
    dY_st_sub = Y_st;
    ddY_st_sub = Y_st;
    [Y_st_sub dY_st_sub ddY_st_sub] = evalNonLin((fit.I_sc(~inds_fixed,:).*fit.A_sc(~inds_fixed,:))*X_ct,fit.f);
    Y_st(inds_fixed,:) = Y_st_fixed;
    Y_st(~inds_fixed,:) = Y_st_sub;
    dY_st(inds_fixed,:) = dY_st_fixed;
    dY_st(~inds_fixed,:) = dY_st_sub;
    ddY_st(inds_fixed,:) = dY_st_fixed;
    ddY_st(~inds_fixed,:) = dY_st_sub;
    %[Y_st dY_st ddY_st] = evalNonLin((fit.I_sc.*fit.A_sc)*X_ct,fit.f);
    % get the current RGC output
    [Z_t dZ_t ddZ_t] = evalNonLin((fit.J_s.*fit.B_s)*Y_st,fit.g);
end


% check for small values, set to 0 for lik and gradient
etol = 10^-6;
iiz = Z_t <= etol;
Z_t(iiz) = etol;
if exist('dZ_t','var'); dZ_t(iiz) = 0; end
if exist('ddZ_t','var'); ddZ_t(iiz) = 0; end

% compute liklihood
loglik = getLogLikSpk(Z_t,R_t);
err = -loglik;

% the gradient has a simpler form when g is exponential, but it doesn't
% depend on the form of f, so we compute that first

if nargout > 1
    switch fit.g.type
        case 'exp'
            % get gradient
            W_s = fit.J_s.*fit.B_s;
            dY_st_W_s = bsxfun(@times,dY_st,W_s');
            grad_A = -(bsxfun(@times,dY_st_W_s,R_t-dZ_t)*X_ct');
            grad_A = grad_A(prinds);
            grad_mu = -sum((Z_t-R_t).*(W_s*dY_st));
            grad = [vector(grad_A); grad_mu];
        otherwise
            Q_t = dZ_t./Z_t;
            W_s = fit.J_s.*fit.B_s;
            dY_st_W_s = bsxfun(@times,dY_st,W_s');
            grad_A = -(bsxfun(@times,dY_st_W_s,(R_t.*Q_t)-dZ_t)*X_ct');
            grad_A = grad_A(prinds);
            if usingSplines
                grad = vector(grad_A);
            else
                grad_mu = -sum((dZ_t-R_t.*(dZ_t./Z_t)).*(W_s*dY_st));
                grad = [vector(grad_A); grad_mu];
            end
    end
end

if nargout > 2
    switch fit.f.type
        case {'rectLinear','rectLinearNeg'} % simpler form when f is rectified
            switch fit.g.type
                case 'exp' % simpler form when g is exponential
                    %% get hessian
                    i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
                    j = 1:size(X_ct,2); j = j(ones(1,1),:);
                    mult1 = X_ct(i,j);
                    mult2 = repmat(dY_st_W_s,[fit.C 1]);
                    all_vec_1 = mult1(prinds,:) .* mult2(prinds,:);
                    H_A = -all_vec_1*bsxfun(@times,all_vec_1,-Z_t)';
                    H_mu = sum((Z_t).*((W_s*dY_st).^2));
                    H_mu_A = -sum(bsxfun(@times,all_vec_1,(Z_t).*(W_s*dY_st)),2);
                    hess = [H_A H_mu_A; H_mu_A', H_mu];
                    
                otherwise % non-exponential g, but f is rectified      
                    % get hessian
                    Q_t = ((Z_t.*ddZ_t)-(dZ_t.*dZ_t))./(Z_t.*Z_t);
                    i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
                    j = 1:size(X_ct,2); j = j(ones(1,1),:);
                    mult1 = X_ct(i,j);
                    mult2 = repmat(dY_st_W_s,[fit.C 1]);
                    all_vec_1 = mult1(prinds,:) .* mult2(prinds,:);
                    H_A = -all_vec_1*bsxfun(@times,all_vec_1,(R_t.*Q_t-ddZ_t))';
                    H_mu = -sum((R_t.*Q_t-ddZ_t).*((W_s*dY_st).^2));
                    H_mu_A = -sum(bsxfun(@times,all_vec_1,(-1*R_t.*Q_t+ddZ_t).*(W_s*dY_st)),2);
                    hess = [H_A H_mu_A; H_mu_A', H_mu];
                    
            end
        otherwise % f is something other than rectified
            %% just going to do this for a generic g, for now
            Q_t = ((Z_t.*ddZ_t)-(dZ_t.*dZ_t))./(Z_t.*Z_t);
            Q_t_2 = dZ_t./Z_t;
            
            i = 1:size(X_ct,1); i = i(ones(1,fit.S),:);
            j = 1:size(X_ct,2); j = j(ones(1,1),:);
            mult1 = X_ct(i,j);
            mult2 = repmat(dY_st_W_s,[fit.C 1]);
            all_vec_1 = mult1(prinds,:) .* mult2(prinds,:);
            %%
            ddY_st_W_s = bsxfun(@times,ddY_st,W_s');
            mult3 = repmat(ddY_st_W_s,[fit.C 1]);
            all_vec_2 = mult1(prinds,:) .* mult3(prinds,:);
            
            H_A_1 = -all_vec_1*bsxfun(@times,all_vec_1,(R_t.*Q_t-ddZ_t))';
            H_A_2 = -all_vec_2*bsxfun(@times,all_vec_2,(R_t.*Q_t_2-dZ_t))';
            H_A = H_A_1 + H_A_2;
            
            if usingSplines
                hess = [H_A];
            else
                % need to fix this...
                H_mu_1 = -sum((R_t.*Q_t-ddZ_t).*((W_s*dY_st).^2));
                H_mu = H_mu_1;
                
                H_mu_A_1 = -sum(bsxfun(@times,all_vec_1,(-1*R_t.*Q_t+ddZ_t).*(W_s*dY_st)),2);
                mult4 = repmat(ddY_st,[fit.C 1]);
                H_mu_A_2 = -sum(bsxfun(@times,mult4(prinds,:),dZ_t-R_t.*Q_t_2),2);
                H_mu_A = H_mu_A_1 + H_mu_A_2;
                
                hess = [H_A H_mu_A; H_mu_A', H_mu];
            end
    end
end

