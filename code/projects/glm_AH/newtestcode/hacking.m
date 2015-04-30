% Dirty code to get up the fact that predicting X with fit Y doesn't work
% 2014-06-01

%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE! 
clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
baseoutput_dir = '/Users/akheitman/NSEM_Home/CrossStim_Performance';

% SETUP cells and experiments, the TYPE of GLM (GLMType) 
BD = NSEM_BaseDirectories;
exptests = [1 2 3 4];
cellselectiontype = 'shortlist';%cellselectiontype = 'debug';
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.specialchange = false;
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
GLMType.CONVEX = true;GLMType.DoubleOpt = true;
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
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

%%
for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    cells
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    SPars = StimulusPars.slv;
        
    
    
    
    
    %%%%  Shorten Block count if using Debug
    if GLMType.debug
        StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:10);
    end
    clear boolean_debug map_type fit_type shead_cellID expname 
    
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = GLMType.map_type; 
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
  
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs; 
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));  
    display(sprintf('Save Directory :  %s', d_save));
    if ~exist(d_save), mkdir(d_save); end
    GLMType.d_save = d_save; 
    
    %% Load Movie and Concatenate the Fitting Section
    clear Main_SolPars Other_SolParss

    [testmovie,inputstats] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
    clear origmatfile
    clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
    %}
 
    %% Load Cell Specific Elements   Spikes and STA
    inputs.exp_nm       = exp_nm; 
    inputs.map_type     = GLMType.map_type; 
    DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs); 
    inputs.stim_type    = GLMType.fit_type
    DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs); 
    
    clear inputs
    
    newdsave = sprintf('%s_newxval',d_save);
    if ~exist(newdsave,'dir'), mkdir(newdsave); end
    
    for i_cell = 1:length(cells)
        clear glm_cellstruct
        
        cid = cells{i_cell};
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);       
        outputdir = sprintf('%s/%s/4sec_%s', baseoutput_dir, GLMType.fitname,exp_nm)  
        if ~exist(outputdir, 'dir'), mkdir(outputdir); end
         
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
        eval(sprintf('load %s/%s.mat', d_save, cell_savename));      
        fittedGLM0 = fittedGLM;
        fittedGLM  = rmfield(fittedGLM,'xvalperformance')

        xvalperformance = eval_xvalperformance_NEW(fittedGLM, SPars, organizedspikes,testmovie,inputstats)
        fittedGLM.xvalperformance = xvalperformance;
        
        printname = sprintf('%s/DiagPlots_%s', newdsave, cell_savename)
        printglmfit(fittedGLM,printname);
        
        eval(sprintf('save %s/%s.mat fittedGLM', newdsave, cell_savename))
        
        
    end
    
end