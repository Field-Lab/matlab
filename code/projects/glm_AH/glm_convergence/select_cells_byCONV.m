%{
 FittedConv.fittedCONV_OFFP{1}
ans = 

       xval_WN_netdecay: 10
     xval_NSEM_netdecay: 10
      train_WN_netdecay: 3.5000
    train_NSEM_netdecay: 1.5000
          xval_WN_t_min: 2.1142
        xval_NSEM_t_min: 4.2283
         train_WN_t_min: 6.0405
       train_NSEM_t_min: 28.1889
         skipminutes_WN: 6.8833
       skipminutes_NSEM: 13.7667
% load
% /Users/akheitman/NSEM_Home/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots/FitConv_Skip25Percent/fittedCONV_1.mat
% /Users/akheitman/NSEM_Home/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots/FitConv_Skip25Percent/fittedCONV_2.mat

% FittedConv.fittedCONV_OFFP{i_cell}.train_WN_netdecay;
% FittedConv.fittedCONV_OFFP{i_cell}.train_NSEM_netdecay;

%}

% AKHEITMAN 2015-03-03
% load /Users/akheitman/NSEM_Home/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots/FitConv_Skip25Percent/fittedCONV_1.mat

clear; close all; clc
load_dir = '/Users/akheitman/NSEM_Home/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots/FitConv_Skip25Percent';

% SET DIRECTORIES
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
raster_scores = allcells;

AA = allcells;  i_exp = 1; i_celltype = 1;
conv_criterion = 2:.5:5;



%%
for i_exp = 1:4
    %%
    clear FittedConv
    eval(sprintf('load %s/fittedCONV_%d.mat', load_dir, i_exp))
    
    AA{i_exp}.conv_criterion = conv_criterion;
    
    OFFP_scores = FittedConv.fittedCONV_OFFP;
    ONP_scores  = FittedConv.fittedCONV_ONP;
    
    ONP_CONV  = NaN(length(allcells{i_exp}.ONP)  , length(conv_criterion));
    OFFP_CONV = NaN(length(allcells{i_exp}.OFFP) , length(conv_criterion));
    
    ONP  = allcells{i_exp}.ONP;
    OFFP = allcells{i_exp}.OFFP;
    
    for i_cell = 1:length(ONP)
        minscore = min(ONP_scores{i_cell}.train_WN_netdecay, ONP_scores{i_cell}.train_WN_netdecay);
        
        good_ind = find(conv_criterion <= minscore);
        bad_ind  = find(conv_criterion > minscore);
        
        ONP_CONV(i_cell,good_ind) = 1;
        ONP_CONV(i_cell,bad_ind)  = 0;
    end
    for i_cell = 1:length(OFFP)
        minscore = min(OFFP_scores{i_cell}.train_WN_netdecay, OFFP_scores{i_cell}.train_WN_netdecay);
        
        good_ind = find(conv_criterion <= minscore);
        bad_ind  = find(conv_criterion > minscore);
        
        OFFP_CONV(i_cell,good_ind) = 1;
        OFFP_CONV(i_cell,bad_ind)  = 0;
    end
    
    AA{i_exp}.ONP_CONV  =  ONP_CONV;
    AA{i_exp}.OFFP_CONV = OFFP_CONV;
        
end

%%

allcells_glmconv = AA;
eval(sprintf('save %s/allcells_glmconv.mat allcells_glmconv', BD.Cell_Selection))