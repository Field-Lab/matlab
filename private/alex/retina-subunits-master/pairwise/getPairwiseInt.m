function out = getPairwiseInt(train,test,ind1,ind2,fit,testDeriv)

% do the fitting
sig1 = train.X_ct(ind1,:)';
sig2 = train.X_ct(ind2,:)';
spks = train.R_t';

%-----------------------------------%
warning('off')
hessfun = @(x,lambda) getOutputs(@() pairwiseNullErrFun(x(:),sig1,sig2,spks,fit),3);
options = optimset('Display','off','GradObj','on','Hessian','user-supplied','HessFcn',hessfun);
%options = optimset('Display','off','GradObj','on');
null_nParams = 3;
null_err_tmp = @(x) pairwiseNullErrFun(x,sig1,sig2,spks,fit);
initPrs = [-0.2 0.2 0.2];
initPrs = rand(size(initPrs))+0.5;
if testDeriv
    DerivCheck_Elts(null_err_tmp,initPrs(:));
    HessCheck(null_err_tmp,initPrs);
    out = [];
end
cons = [-10 zeros(1,null_nParams-1)];
null_params = fmincon(null_err_tmp,initPrs,[],[],[],[],cons,[],[],options);

%-----------------------------------%
hessfun = @(x,lambda) getOutputs(@() pairwiseAltErrFun(x(:),sig1,sig2,spks,fit),3);
options = optimset('Display','off','GradObj','on','Hessian','user-supplied','HessFcn',hessfun);
%options = optimset('Display','off','GradObj','on');
alt_nParams = 5;
alt_err_tmp = @(x) pairwiseAltErrFun(x,sig1,sig2,spks,fit);
initPrs = [-0.2 0.2 0.2 0.1 0.1];
initPrs = rand(size(initPrs))+0.5;
if testDeriv
    DerivCheck_Elts(alt_err_tmp,initPrs(:));
    HessCheck(alt_err_tmp,initPrs);
    out = [];
    return
end
cons = [-10 zeros(1,alt_nParams-1)];
alt_params = fmincon(alt_err_tmp,initPrs,[],[],[],[],cons,[],[],options);


%-----------------------------------%
% get errors on training set
null_err = null_err_tmp(null_params);
alt_err = alt_err_tmp(alt_params);

% get errors on testing set
sig1_test = test.X_ct(ind1,:)';
sig2_test = test.X_ct(ind2,:)';
spks_test = test.R_t';
null_err_tmp_test = @(x) pairwiseNullErrFun(x,sig1_test,sig2_test,spks_test,fit);
alt_err_tmp_test = @(x) pairwiseAltErrFun(x,sig1_test,sig2_test,spks_test,fit);
null_err_test = null_err_tmp_test(null_params);
alt_err_test = alt_err_tmp_test(alt_params);

% collect output
out.null_params = null_params;
out.null_err = -null_err;
out.null_err_test = -null_err_test;
out.alt_params = alt_params;
out.alt_err = -alt_err;
out.alt_err_test = -alt_err_test;

% get test statistics
out.testStat = out.alt_err - out.null_err;
out.testStat_test = out.alt_err_test - out.null_err_test;

warning('on')