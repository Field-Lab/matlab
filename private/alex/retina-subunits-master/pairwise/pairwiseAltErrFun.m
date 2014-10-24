function [nll grad hess] = pairwiseAltErrFun(prs,sig1,sig2,r,fit)

constants = prs(1);
weights = prs(2:5);

[f1 df1] = evalNonLin(weights(1)*sig1,fit.f);
[f2 df2] = evalNonLin(weights(2)*sig2,fit.f);
[f3 df3] = evalNonLin(weights(3)*sig1 + weights(4)*sig2,fit.f);

b = (constants(1));

if strcmp(fit.g.type,'exp')
    [g] = evalNonLin(f1 + f2 + f3 + b,fit.g);
else
    [g dg] = evalNonLin(f1 + f2 + f3 + b,fit.g);
end

etol = 10^-6;
iiz = g <= etol;
g(iiz) = etol;

if ~strcmp(fit.g.type,'exp')
    dg(iiz) = 0;
end

loglik = getLogLikSpk(g,r);
nll = -loglik;

if nargout > 1
    if strcmp(fit.g.type,'exp')
        grad1 = -(sum(r.*df1.*sig1) - sum(g.*df1.*sig1));
        grad2 = -(sum(r.*df2.*sig2) - sum(g.*df2.*sig2));
        grad3 = -(sum(r.*df3.*sig1) - sum(g.*df3.*sig1));
        grad4 = -(sum(r.*df3.*sig2) - sum(g.*df3.*sig2));
        gradb = -(sum(r) - sum(g));
        grad = [gradb grad1 grad2 grad3 grad4];
    else
        q = r.*(dg./g);
        grad1 = -(sum(q.*df1.*sig1) - sum(dg.*df1.*sig1));
        grad2 = -(sum(q.*df2.*sig2) - sum(dg.*df2.*sig2));
        grad3 = -(sum(q.*df3.*sig1) - sum(dg.*df3.*sig1));
        grad4 = -(sum(q.*df3.*sig2) - sum(dg.*df3.*sig2));
        gradb = -(sum(q) - sum(dg));
        grad = [gradb grad1 grad2 grad3 grad4];
    end
end

% this hessian is only correct for rectified linear and exponential
if nargout > 2
    tmp1 = sig1.*df1;
    tmp2 = sig2.*df2;
    tmp13 = sig1.*df3;
    tmp23 = sig2.*df3;
    
    H11 = sum(g.*(tmp1).^2);
    H22 = sum(g.*(tmp2).^2);
    H33 = sum(g.*(tmp13).^2);
    H44 = sum(g.*(tmp23).^2);
    
    H12 = sum(g.*tmp1.*tmp2);
    H13 = sum(g.*tmp1.*tmp13);
    H14 = sum(g.*tmp1.*tmp23);
    H23 = sum(g.*tmp2.*tmp13);
    H24 = sum(g.*tmp2.*tmp23);
    H34 = sum(g.*tmp13.*tmp23);
    
    H1b = sum(g.*tmp1);
    H2b = sum(g.*tmp2);
    H3b = sum(g.*tmp13);
    H4b = sum(g.*tmp23);
    Hbb = sum(g);
    hess = [Hbb H1b H2b H3b H4b; 
            H1b H11 H12 H13 H14;
            H2b H12 H22 H23 H24;
            H3b H13 H23 H33 H34;
            H4b H14 H24 H34 H44];
end