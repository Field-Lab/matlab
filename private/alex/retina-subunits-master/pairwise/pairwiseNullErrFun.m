function [nll grad hess] = pairwiseNullErrFun(prs,sig1,sig2,r,fit)

constants = prs(1);
weights = prs(2:3);

[f1 df1] = evalNonLin(weights(1)*sig1,fit.f);
[f2 df2] = evalNonLin(weights(2)*sig2,fit.f);

b = (constants(1));

if strcmp(fit.g.type,'exp');
    [g ] = evalNonLin(f1 + f2 + b,fit.g);
else
    [g dg] = evalNonLin(f1 + f2 + b,fit.g);
end

etol = 10^-6;
iiz = g <= etol;
g(iiz) = etol;

if ~strcmp(fit.g.type,'exp');
    dg(iiz) = 0;
end

loglik = getLogLikSpk(g,r);
nll = -loglik;

if nargout > 1
    if strcmp(fit.g.type,'exp')
        grad1 = -(sum(r.*df1.*sig1) - sum(g.*df1.*sig1));
        grad2 = -(sum(r.*df2.*sig2) - sum(g.*df2.*sig2));
        gradb = -(sum(r) - sum(g));
        grad = [gradb grad1 grad2];
    else
        q = r.*(dg./g);
        grad1 = -(sum(q.*df1.*sig1) - sum(dg.*df1.*sig1));
        grad2 = -(sum(q.*df2.*sig2) - sum(dg.*df2.*sig2));
        gradb = -(sum(q) - sum(dg));
        grad = [gradb grad1 grad2];
    end
end

% this hessian is only correct for rectified linear and exponential
if nargout > 2
    tmp1 = sig1.*df1;
    tmp2 = sig2.*df2;
    H12 = sum(g.*tmp2.*tmp1);
    H11 = sum(g.*(tmp1).^2);
    H22 = sum(g.*(tmp2).^2);
    H1b = sum(g.*tmp1);
    H2b = sum(g.*tmp2);
    Hbb = sum(g);
    hess = [Hbb H1b H2b; 
            H1b H11 H12;
            H2b H12 H22];
end
    