clear; close all; clc
basefolder = '/Users/akheitman/NSEM_Home/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL';

for i_exp = exps
    for i_stimtype = stimtypes 
        for i_celltype = celltypes


cid = 5086;
cell_savename = 'OFFPar_5086'
exp_nm = '2012-08-09-3';
pct_change = [5:10:95];pct_change = [pct_change, 100];
i_pct = 1;


%%

% Identify Time Scales of Fitting
conv.WN.fit_seconds = NaN(1,length(pct_change));
conv.NSEM.fit_seconds = NaN(1,length(pct_change));
for i_pct = 1:length(pct_change)
    pct = pct_change(i_pct);
    [StimulusPars] = Directories_Params_v23_Conv_hack(exp_nm,'WN',pct/100);
    slv_WN = StimulusPars.slv;
    conv.WN.fit_seconds(i_pct) = (length(slv_WN.FitBlocks) * length(slv_WN.fitframes)) / 120;
    
    [StimulusPars] = Directories_Params_v23_Conv_hack(exp_nm,'NSEM',pct/100);
    slv_NSEM = StimulusPars.slv;
    conv.NSEM.fit_seconds(i_pct) = (length(slv_NSEM.FitBlocks) * length(slv_NSEM.fitframes)) / 120;
end
conv.WN.fit_minutes   = conv.WN.fit_seconds/60;
conv.NSEM.fit_minutes = conv.NSEM.fit_seconds/60;
conv.WN.xval_bps              = zeros(1, length(pct_change));
conv.WN.objective_val         = zeros(1, length(pct_change)); 
conv.NSEM.xval_bps            = zeros(1, length(pct_change));
conv.NSEM.objective_val       = zeros(1, length(pct_change));

conv.WN.fit_p           = cell(length(pct_change),1 )
conv.NSEM.fit_p         = cell(length(pct_change),1 )
conv.WN.linearfilters   = cell(length(pct_change),1 );
conv.NSEM.linearfilters = cell(length(pct_change),1 );

% Roughly One Second Per Cell   LOAD WN AND NSEM VALS %
for i_pct = 1:length(pct_change)
    pct = pct_change(i_pct);
    
    WN_dir   = sprintf('ChangeParams_Fit_Convergence_%dPct/WN_mapPRJ/%s', pct, exp_nm);
    NSEM_dir = sprintf('ChangeParams_Fit_Convergence_%dPct/NSEM_mapPRJ/%s', pct, exp_nm);
    
    for i_fit = 1:2
        if i_fit == 1, eval(sprintf('load %s/%s/%s.mat', basefolder, WN_dir, cell_savename));   Z = conv.WN; end
        if i_fit == 2, eval(sprintf('load %s/%s/%s.mat', basefolder, NSEM_dir, cell_savename)); Z =conv.NSEM; end
        
        Z.xval_bps(i_pct)       = fittedGLM.xvalperformance.logprob_glm_bpspike;
        Z.objective_val(i_pct)  = fittedGLM.rawfit.objective_val;
        
        if i_fit == 1, conv.WN   = Z; end
        if i_fit == 2, conv.NSEM = Z; end
    end
end
%%
%{
figure;
panels = 5
for i_fit = 1:2
    if i_fit == 1, Z = conv.WN; fit = 'WN';  string_color = 'b';   end
    if i_fit == 2, Z = conv.NSEM; fit = 'NSEM'; string_color = 'r'; end
    
    for i_panel = 1:panels
        switch i_panel
            case 1
                measure = 'XVAL BPS'; xvals = Z.fit_minutes; yvals = Z.xval_bps;
            case 2
                measure = 'ObjVal'; xvals = Z.fit_minutes; yvals = Z.objective_val;
            case 3
                measure = 'ObjVal/fittime'; xvals = Z.fit_minutes; yvals = Z.objective_val./Z.fit_minutes;
            case 4
                measure = 'Increment ObjVal'; xvals = Z.fit_minutes(2:end); yvals = diff(Z.objective_val);
            case 5
                measure = 'Deriv ObjVal'; 
                xvals = Z.fit_minutes(2:end-1);
                center_point = Z.objective_val(2:end-1);
                before_point = Z.objective_val(1:end-2); after_point = Z.objective_val(3:end);
                backward_deriv = (center_point-before_point)./ diff(Z.fit_minutes(1:end-1));
                forward_deriv = (after_point-center_point)./ diff(Z.fit_minutes(2:end));
                yvals = .5 * backward_deriv + .5*forward_deriv;
        end
        
       subplot(2,panels, (i_panel+(i_fit-1)*panels) ); hold on; 
        set(gca,'fontsize',8); title(sprintf('%s: %s', fit,measure));
        plot(xvals,yvals,string_color);
        plot(xvals,yvals,'k.'); hold off
    end
end
%}        
        

%%
        
        
    
    
    
    

