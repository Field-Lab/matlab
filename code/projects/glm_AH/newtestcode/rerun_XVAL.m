% Dirty code to get up the fact that predicting X with fit Y doesn't work
% 2014-06-01

%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE! 
clear; close all; clear all; clc
% Use this for the convergence check!


% SETUP cells and experiments, the TYPE of GLM (GLMType) 
BD = NSEM_BaseDirectories;
exptests = [3];

%cellselectiontype = 'shortlist';%cellselectiontype = 'debug';
cellselectiontype = 'all';
clear GLMType
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.nullpoint = 'mean'; 
GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.fit_type = 'WN';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = true
%

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
%%
for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    SPars = StimulusPars.slv;
    
    if strcmp(cellselectiontype, 'all')
        [exp_nm,cells,expname,badcells]  = cell_list( expnumber, cellselectiontype);
        clear cells
        [~, ~, datarun_slv_WN datarun_mas]   = Directories_Params_v23_Conv(exp_nm, 'WN', 'mapPRJ',1);
        [~, ~, datarun_slv_NSEM datarun_mas] = Directories_Params_v23_Conv(exp_nm, 'NSEM', 'mapPRJ',1);
        ONP  = intersect(datarun_slv_WN.cell_types{1}.cell_ids , datarun_slv_NSEM.cell_types{1}.cell_ids);
        OFFP = intersect(datarun_slv_WN.cell_types{2}.cell_ids , datarun_slv_NSEM.cell_types{2}.cell_ids);
        for i_cell = 1:length(badcells);
            ONP(find(ONP == badcells(i_cell) )) = [];
            OFFP(find(OFFP == badcells(i_cell) )) = [];
        end
        %cells_vec = [ONP OFFP];
        cells_vec = [ONP];
        cells = cell(1,length(cells_vec));
        for i_cell = 1:length(cells);
            cells{i_cell} = cells_vec(i_cell);
        end
    end
    
    
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
    %%% Load Stimulus   -- insert more frame cutting here!    
    %[blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
    
    [testmovie, inputstats]   = loadmoviematfile(exp_nm, GLMType.fit_type, GLMType.cone_model,'testmovie');
    clear origmatfile
    clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
    %}
 
    %% Load Cell Specific Elements   Spikes and STA
    inputs.exp_nm       = exp_nm; 
    inputs.map_type     = GLMType.map_type; 
    DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs); 

    inputs.stim_type    = GLMType.fit_type;
    DirPars.organizedspikesdir= NSEM_secondaryDirectories('organizedspikes_dir', inputs); 

    

    for i_cell = 1:length(cells)
        clear glm_cellstruct
        
        cid = cells{i_cell};
        display(sprintf('Working on cell %d', cid));
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);       
        eval(sprintf('load %s/%s.mat', d_save, cell_savename));
       
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
        fittedGLM = rmfield(fittedGLM,'xvalperformance');        
     
        xvalperformance = eval_xvalperformance_WithSMOOTHING(fittedGLM, SPars, organizedspikes,testmovie, inputstats,[1 2 4 6 8]);
        
        fittedGLM.xvalperformance = xvalperformance;
        eval(sprintf('save %s/%s_generalXVAL.mat xvalperformance', d_save, fittedGLM.cellinfo.cell_savename));
        %printname = sprintf('%s/DiagPlots_%s', d_save,fittedGLM.cellinfo.cell_savename);
        %printglmfit(fittedGLM,printname)
        
    end
    

end
