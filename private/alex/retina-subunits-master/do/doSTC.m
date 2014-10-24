function [r2_LN r2_STC] = doSTC(celldat,train,test)

% given a celldat structure, training data, and test data,
% do an STC analysis and save out a figure with the
% STA and STC filters

% get R2 for an LN fit
fit_LN = getSTC(train,0);
fit_LN.g.type = 'logexp1';
fit_LN.model = 'ln';
fit_LN.error = 'loglik';
fit_LN.g.p = [0 1];
fit_LN.D = 1;
fit_LN.d_s = [fit_LN.k'];
fit_LN = fitLN(train,fit_LN,0);
fit_LN.g = initSpline(fit_LN.g);
fit_LN = fitG(train,fit_LN);
out_LN = evalFit(test,fit_LN);

% get R2 for an STC fit
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

% save the LN file
res = [];
res.fit_LN = fit_LN;
res.out_LN = out_LN;
if ~isdir(strcat(celldat.analPath,'/LN/'))
  mkdir(strcat(celldat.analPath,'/LN/'));
end
saveName = sprintf('%g-LN',celldat.rgcId);
save(fullfile(celldat.analPath,'LN',strcat(saveName,'-fit.mat')),'res');

% save the STC file
res = [];
res.fit_STC = fit_STC;
res.out_STC = out_STC;
if ~isdir(strcat(celldat.analPath,'/STC/'))
  mkdir(strcat(celldat.analPath,'/STC/'));
end
saveName = sprintf('%g-STC',celldat.rgcId);
save(fullfile(celldat.analPath,'STC',strcat(saveName,'-fit.mat')),'res');

