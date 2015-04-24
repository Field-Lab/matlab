% Dirty Code to get cross validated bits per spike scores up 
% AKHeitman 2014-11-23


%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE! 
clear; close all; clc
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire

t_bin        = .00083275;
%type         = 'ViktorSpike'; 
binscales    = [1 2 5 10 25 50 100];
%binscales    = [1];
% SETUP cells and experiments, the TYPE of GLM (GLMType) 
BD = NSEM_BaseDirectories;
%exptests = [1 3];
exptests = [2 4];
selectiontype = 'shortlist';%cellselectiontype = 'debug';
cellselectiontype = 'all';
clear GLMType
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.nullpoint = 'mean'; 
GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = true


GLMType.specialchange = false;
%}
%{
GLMType.specialchange = true;
GLMType.specialchange_name = 'Fit_Convergence';
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.DoubleOpt = true;
GLMType.DoubleOpt_Manual = true;
%}

GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.func_sname = 'glmwrap_23';
GLMType.fullmfilename =mfilename('fullpath'); 
i_exp = 1; i_cell = 1;
GLMType.fitname  = GLM_fitname(GLMType);
if GLMType.specialchange && strcmp(GLMType.specialchange_name,'Fit_Convergence')
    GLMType.fitname = sprintf('%s_100Pct', GLMType.fitname);
end

troubleshoot.doit    = true;
troubleshoot.plotdir = '/Users/akheitman/Matlab_code/troubleshooting_plots'
%troubleshoot.plotdir = BD.GLM_troubleshootplots  % temporarily change to
troubleshoot.name    = 'singleopt';
i_exp = 1; i_celltype = 1; i_cell = 1; i_time = 1;
%%

for i_stim = 1:2
    if i_stim == 1, GLMType.fit_type = 'WN'; end
    if i_stim == 2, GLMType.fit_type = 'NSEM'; end
    GLMType.fitname  = GLM_fitname(GLMType)
for i_exp = exptests
    %% 
    exp_nm = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    ONP = allcells{i_exp}.ONP;
    OFFP = allcells{i_exp}.OFFP;
    cells = [ONP , OFFP];
    %%%%  Shorten Block count if using Debug
    if GLMType.debug
        StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:10);
    end
    clear boolean_debug map_type fit_type shead_cellID expname 
    
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
    % Get directories taken care of %     
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = GLMType.map_type; 
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
    inputs.stim_type    = GLMType.fit_type;
  
    dir.fittedGLM = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs; 
    dir.rastCompare = sprintf('%s/rastCompare', dir.fittedGLM);
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));  
    display(sprintf('Save Directory :  %s', dir.rastCompare));
    if ~exist(dir.rastCompare), mkdir(dir.rastCompare); end
   
    %%
    for i_celltype = [1,2]
        if i_celltype == 1, ctype = 'ONPar'; cells = ONP; end
        if i_celltype == 2, ctype = 'OFFPar'; cells = OFFP; end
        
        for i_cell = 1:length(cells)
            %%
            clear rastCompare_error
            cid = cells(i_cell);
            display(sprintf('Working on %s : %d', ctype, cid));
            cell_savename = sprintf('%s_%d', ctype, cid);    
            eval(sprintf('load %s/%s.mat', dir.fittedGLM, cell_savename));
            
            rastRec       = fittedGLM.xvalperformance.rasters.recorded; 
            rastSim       = fittedGLM.xvalperformance.rasters.glm_sim;
            rastCompare_error.cid = cid;
            rastCompare_error.cell_savename = cell_savename;
            rastCompare_error.exp_nm = exp_nm;
            rastCompare_error.ViktorDistperSpike.bytime    = cell(1,length(binscales));
            rastCompare_error.ViktorDistperSpike.timescale = 2* t_bin * binscales; % 2 to bring close to register with 
            
            rastCompare_error.L2_error.bytime     = cell(1,length(binscales));
            rastCompare_error.L2_error.smoothbins = binscales;
            rastCompare_error.L2_error.bintime    = t_bin;
            rastCompare_error.L2_error.note = 'Convolve spikes by normal distrbution with sigma at time scale of smoothbins to get rate';
            for i_time = 1:length(binscales)
                % Viktor Spikes
                t_cost = 2 * t_bin * binscales(i_time);
                [vksp_mean , vksp_std] = eval_xvalperformance_ViktorDistperSpike(rastRec,rastSim,t_bin,t_cost);
                rastCompare_error.ViktorDistperSpike.bytime{i_time}.mean = vksp_mean;
                rastCompare_error.ViktorDistperSpike.bytime{i_time}.std = vksp_std;
                 % L2 vector norms
                smoothbins = binscales(i_time);
                L2_Error = eval_xvalperformance_L2error(rastRec,rastSim, t_bin, smoothbins);
                rastCompare_error.L2_error.bytime{i_time}.L2_Error = L2_Error;
            end
            eval(sprintf('save %s/%s_RastComparisonMetrics.mat rastCompare_error', dir.rastCompare, fittedGLM.cellinfo.cell_savename));

        end
    end

end

end
