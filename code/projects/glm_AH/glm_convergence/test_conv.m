% Trying Convergence Cut
% AKHEITMAN 2015-02-08

% Scale Convergence to [-0, -1];
% Then choose differe

%%
clear; clc;
BD = NSEM_BaseDirectories;
basefolder   = sprintf('%s/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots', BD.NSEM_home);
celltypes = [1 2];
exps = [1 2 3 4];
st_index = 3; % 25th percentile of Data
i_exp = 1;
i_celltypes = 1;


for i_exp = exps
    eval(sprintf('load %s/expTrain_Conv_%d_shortlist.mat',basefolder,i_exp))
    TC = expTrain_Conv;
    
    % SCALE TIME TO 1
    WN.timing.delta_t =  TC.NSEM_fit_minutes(end) - TC.NSEM_fit_minutes(st_index);
    WN.timing.lags    = (TC.NSEM_fit_minutes(st_index:end) - TC.NSEM_fit_minutes(st_index)) / WN.timing.delta_t;
    NSEM.timing.delta_t =  TC.WN_fit_minutes(end) - TC.WN_fit_minutes(st_index);
    NSEM.timing.lags    = (TC.WN_fit_minutes(st_index:end) - TC.WN_fit_minutes(st_index)) / NSEM.timing.delta_t;

    % ASSUME EXPONENTIAL(-1/tau)  FIND CORRESPONDING TOP POINT
    taus_unitless = 1:10;
    startpoint    = 1./(1-exp(-1./taus_unitless));
    
    for i_celltype = celltypes
        if i_celltype == 1; cellgroup = TC.ONP;  celltype = 'ONPar'; train_data = TC.test_train_ONP; end
        if i_celltype == 2; cellgroup = TC.OFFP; celltype = 'OFFPar'; train_data = TC.test_train_OFFP;end
        test_train.WN    = cell(length(cellgroup),1);
        test_train.NSEM  = cell(length(cellgroup),1); 
        for i_cell = 1:length(cellgroup)
            
            
            
            NSEM.raw_vals          = TC.test_train_OFFP.NSEM{i_cell}.objval;
            normedfit_fromend = (raw_vals(st_index:end) - raw_vals(st_index)) / (raw_vals(st_index) - raw_vals(end));
        end
    end
end



%%
close all
l2error = zeros(1,length(taus))
for i_tau = 1:length(taus)
    conv_vals = normedfit_fromend + multiplefront(i_tau);
    fit_vals = (multiplefront(i_tau)) * exp(-lags / taus(i_tau));
    
    figure; hold on; plot(conv_vals,'r'); plot(fit_vals,'g');
    l2error(i_tau) = norm(conv_vals - fit_vals');
end





%%






%%



