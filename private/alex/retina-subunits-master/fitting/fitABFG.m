function [fit_SUB negloglik r2] = fitABFG(dat,train,test,structure,fitOld)

fit_SUB.error = 'loglik';
fit_SUB.model = 'subunit';
fit_SUB.constraints.B_s = 'none';
fit_SUB.constraints.A_sc = 'L1';
fit_SUB.C = dat.nConesFit;
fit_SUB.locs_c = dat.locs_c;

fit_SUB.I_sc = structure;

fit_SUB.S = size(fit_SUB.I_sc,1);
fit_SUB.J_s = ones(1,size(fit_SUB.I_sc,1));
fit_SUB.g.type = 'logexp1';
fit_SUB.g.p = [0 1];
fit_SUB.B_s = normMat(ones(1,fit_SUB.S),'L1',2);
fit_SUB.A_sc = normMat(ones(fit_SUB.S,fit_SUB.C).*fit_SUB.I_sc,'L1',2);

if exist('fitOld','var') && ~isempty(fitOld)
    fit_SUB.f = fitOld.f;
    fit_SUB = fitF(train,fit_SUB,1);
else
    if dat.polarity == 1
    fit_SUB.f.type = 'rectLinear';
    end
    if dat.polarity == -1
    fit_SUB.f.type = 'rectLinearNeg';
    end
    fit_SUB.f.p = 0;
    fit_SUB.f = initSpline(fit_SUB.f);
    fit_SUB = fitF(train,fit_SUB,1);
end
    
%%
nIter = 2;
for iIter = 1:nIter 
    fit_SUB = fitB(train,fit_SUB,0);
    if fit_SUB.S ~= fit_SUB.C || sum(vector(fit_SUB.I_sc ~= eye(fit_SUB.C)))
        fit_SUB = fitA(train,fit_SUB,0);
    end
    fit_SUB = fitF(train,fit_SUB);
end
%%
fit_SUB.g = initSpline(fit_SUB.g);
fit_SUB = fitG(train,fit_SUB);

out_SUB = evalFit(train,fit_SUB);

negloglik = out_SUB.negloglik;
r2 = out_SUB.r2;