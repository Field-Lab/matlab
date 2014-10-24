function doSubunit(celldat,train,test,mode,init_orig);

%%
fit_SUB.error = 'loglik';
fit_SUB.model = 'subunit';
fit_SUB.constraints.B_s = 'none';
fit_SUB.constraints.A_sc = 'L1';
fit_SUB.C = celldat.nConesFit;
fit_SUB.locs_c = celldat.locs_c;

switch mode
    case 'singletons'
        init = mkSubAssgn(fit_SUB.C,'singletons');
        saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode);
    case {'greedy','greedy-local'}
        init = init_orig;
        saveName = sprintf('%g-subunit-%s',celldat.rgcId,mode); 
end
fit_SUB.I_sc = init;
saveName = strcat(saveName,'-',num2str(celldat.percent));

if isfield(celldat,'simopts')
    saveName = strcat(saveName,'-',sprintf('d%g',celldat.simopts.train_dur),'-',num2str(celldat.simopts.iter),'-',celldat.simopts.mode);
end

fprintf('(doSubunit) saving %s\n',saveName);

%%
fit_SUB.S = size(fit_SUB.I_sc,1);
fit_SUB.J_s = ones(1,size(fit_SUB.I_sc,1));
fit_SUB.g.type = 'logexp1';
fit_SUB.g.p = [0 1];
fit_SUB.B_s = normMat(ones(1,fit_SUB.S),'L1',2);
fit_SUB.A_sc = normMat(ones(fit_SUB.S,fit_SUB.C).*fit_SUB.I_sc,'L1',2);

if celldat.polarity == 1
fit_SUB.f.type = 'rectLinear';
end
if celldat.polarity == -1
fit_SUB.f.type = 'rectLinearNeg';
end
fit_SUB.f.p = 0;
fit_SUB.f = initSpline(fit_SUB.f);
fit_SUB = fitF(train,fit_SUB,1);
nIter = 3;
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
%%
out_SUB = evalFit(test,fit_SUB);
%%
fit_LN = getSTC(train,0);
fit_LN.model = 'ln';
fit_LN.error = 'loglik';
fit_LN.g.type = 'logexp1';
fit_LN.g.p = [0 1];
fit_LN.d_s = [fit_LN.k'];
fit_LN = fitLN(train,fit_LN,0);
fit_LN.g = initSpline(fit_LN.g);
fit_LN = fitG(train,fit_LN);
out_LN = evalFit(test,fit_LN);
%%
fit_STC = getSTC(train,0);
fit_STC.model = 'stc';
fit_STC.error = 'loglik';
fit_STC.g.type = 'logexp1';
fit_STC.g.p = [0 1];
fit_STC.D = size(fit_STC.V,2);
fit_STC.d_s = [0.6 0.1*ones(1,fit_STC.D-1)];
fit_STC = fitSTC(train,fit_STC,0);
fit_STC.g = initSpline(fit_STC.g);
fit_STC = fitG(train,fit_STC);
out_STC = evalFit(test,fit_STC);

%%
res.fit_SUB = fit_SUB;
res.out_SUB = out_SUB;
res.fit_LN = fit_LN;
res.out_LN = out_LN;
res.fit_STC = fit_STC;
res.out_STC = out_STC;

if ~isdir(strcat(celldat.analPath,'/subunit/'))
    mkdir(strcat(celldat.analPath,'/subunit/'));
end
save(fullfile(celldat.analPath,'subunit',strcat(saveName,'-fit.mat')),'res');
