function [simdat simtrain simtest] = getSimulation(celldat,train,test,fit,train_dur,test_dur)

simdat = celldat;
simdat.simopts.train_dur = train_dur;
simdat.simopts.test_dur = test_dur;
simdat.analPath = strcat(celldat.mainAnalPath,'simulation-',celldat.dataSet);
simdat.figurePath = strcat(celldat.mainFigurePath,'simulation-',celldat.dataSet);

% get the stim
simtrain.X_ct = train.X_ct(:,1:simdat.simopts.train_dur);
simtest.X_ct = test.X_ct(:,1:simdat.simopts.test_dur);
% get the subunit responses
[simtrain.Y_st] = evalNonLin((fit.I_sc.*fit.A_sc)*simtrain.X_ct,fit.f);
[simtest.Y_st] = evalNonLin((fit.I_sc.*fit.A_sc)*simtest.X_ct,fit.f);
% get the current RGC output
[simtrain.Z_t] = evalNonLin((fit.J_s.*fit.B_s)*simtrain.Y_st,fit.g);
[simtest.Z_t] = evalNonLin((fit.J_s.*fit.B_s)*simtest.Y_st,fit.g);

% simulate spikes, clipping small values
simtrain.R_t = poissrnd(max(simtrain.Z_t,0));
simtest.R_t = poissrnd(max(simtest.Z_t,0));
