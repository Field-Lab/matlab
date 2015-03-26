% Will modernize to be flexbile in how we call cells
% Version 0 worked up till glm_AH_5_3
% 

% Wrap no number done on 2015-01-06
% Wrap_func should be synced with Wrap 2015-01-11
% Wrap 25 started 2014-12-13
% Should just be cleaner all around

% Really is all the ugly loading/cellpicking etc.


% AKHeitman 2014-03-27
% AKHeitman 2014-04-07  -emphasize role as parameter independet loading
%                       -setting directories
%                       -dicatating parameters
%                       -no computations
% Calls:
% GLM_fitname
% cell_list
% Directories_Params_v23
% NSEM_secondaryDirectories
% loadmoviematfile
% concat_fitmovie_fromblockedcell
% findcelltype
% concat_fitspikes_fromorganizedspikes
% visionSTA_to_xymviCoord
%%% 
% glm_execute
%%%%

% Fit Dependent
% Parameter Independent
% Organize fit movie and spike times, directories etc.
% DO ALL LOADING AND SAVING HERE!!!
% All stimulus parameters and blocks everything should get dealt with here


% movie, spikes, block structure, GLMType
% the Convtest calling sequence needs to get worked out.. but otherwise ok


%%%% BELOW ARE SOME EXAMPLE COMMAND SEQUENCES %%%%%
%{
clear
exps = [1]; 
stimtypes = [1];
celltypes = [1];
changes_cell{1}.type = 'debug';
changes_cell{1}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)

clear
exps = [4]; 
stimtypes = [2];
celltypes = [1];
%changes_cell{1}.type = 'debug';
%changes_cell{1}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes)


clear
exps = [1 2 3]; 
stimtypes = [1];
celltypes = [1];
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'rk2';
%changes_cell{2}.type = 'cone_model';
%changes_cell{2}.name = 'rieke_linear';
%changes_cell{2}.type = 'debug';
%changes_cell{2}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
exps = [4]; 
stimtypes = [2];
celltypes = [1];
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'fixedSP-ConductanceBased';
%changes_cell{2}.type = 'debug';
%changes_cell{2}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)

clear
exps = [3]; 
stimtypes = [1];
celltypes = [1];
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'rk2';
%changes_cell{2}.type = 'debug';
%changes_cell{2}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)


exps = [1]; 
stimtypes = [1 2];
celltypes = [2];
changes_cell{1}.type = 'debug';
changes_cell{1}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)

changes_cell{1}.type = 'debug';
changes_cell{1}.name = 'true';
glmwrap_func(1,1,1,changes_cell)


exps = [3 4];
celltypes = [1 2];
stimtypes = 2
changes_cell{1}.type = 'PostSpikeFilter';
changes_cell{1}.name = 'OFF';
%changes_cell{2}.type = 'debug';
%changes_cell{2}.name = 'true';
glmwrap_func(exps,stimtypes,celltypes,changes_cell)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%changes_cell{1}.type = 'filter_mode';
%changes_cell{1}.name = 'rk2-ConductanceBased';
%changes_cell{1}.name = 'rk2-ConductanceBased';
%function glmwrap(exps,stimtypes,celltypes,changes_cell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE! 

function glmwrap_func(exps,stimtypes,celltypes,changes_cell)

BD = NSEM_BaseDirectories;
runoverexisting  = true;
shortlist = false;


reverseorder = false;


if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end

GLMType.fitname    = GLM_fitname(GLMType); 
GLMType.func_sname = 'glmwrap';
GLMType.fullmfilename =mfilename('fullpath'); 
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
i_exp = 1; i_cell = 1;
troubleshoot.doit    = true;


%%
for i_exp = exps
    for i_stimtype = stimtypes 
        for i_celltype = celltypes
            %% Cells / Directories / Parameters
            
            % DEFINE GLMTYPE 
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            GLMType.fit_type = stimtype;
             
            
            % LOAD CELLS AND EXP DEPENDENT VARIABLES
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            [StimulusPars Dirs datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
            if GLMType.debug
                StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:2);
            end

    
            % DIRECTORIES  
            secondDir.exp_nm    = exp_nm; 
            secondDir.map_type  = GLMType.map_type; 
            secondDir.stim_type = GLMType.fit_type;
            secondDir.fitname   = GLMType.fitname;
            Dirs.fittedGLM_savedir  = NSEM_secondaryDirectories('savedir_GLMfit', secondDir);
            Dirs.WN_STAdir          = NSEM_secondaryDirectories('WN_STA', secondDir); 
            Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
            if ~exist(Dirs.fittedGLM_savedir), mkdir(Dirs.fittedGLM_savedir); end         
            display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));  
            display(sprintf('Save Directory :  %s', Dirs.fittedGLM_savedir));
            
            %% Stimulus

            % LOAD BLOCKED STIMULUS IN UINT8 FORM 
            [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
            [testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
            GLMType.fitmoviefile = origmatfile;
            concat_fitmovie      = concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv);
            clear origmatfile
            
            % Hack!! get rid of it 2014-12-15%
            if shortlist,  [~,cellgroup,~]  = cell_list(i_exp, 'shortlist'); cellgroup = cell2mat(cellgroup);  end
            
            
            if exist('reverseorder','var') && reverseorder, cellgroup = fliplr(cellgroup); end
            
            %% Fit GLM for each cell
            for i_cell = 1:length(cellgroup)
                clear glm_cellinfo fittedGLM cell_savename
                cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
                
                % Hack!! get rid of it 2014-12-15%
                if shortlist,  [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);  end
                
                
                if ~exist(sprintf('%s/%s.mat', Dirs.fittedGLM_savedir,cell_savename)) || (exist('runoverexisting','var') && ~runoverexisting)
                    % INITIALIZE GLM_CELLINFO
                    glm_cellinfo.cid           = cid;
                    glm_cellinfo.exp_nm        = exp_nm;
                    glm_cellinfo.celltype      = celltype;
                    glm_cellinfo.cell_savename = cell_savename;
                    glm_cellinfo.fitname       = GLMType.fitname;
                    glm_cellinfo.computedtstim = StimulusPars.slv.computedtstim;
                
                    % CONCATENATE SPIKES
                    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                    spikesconcat.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                    
                    % LOAD STA FROM MASTER DATASET
                    eval(sprintf('load %s/STAandROI_%s.mat STAandROI', Dirs.WN_STAdir, cell_savename));
                    master_idx         = find(datarun_mas.cell_ids == cid);
                    stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
                    stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
                    slvdim.height      = StimulusPars.slv.height; 
                    slvdim.width       = StimulusPars.slv.width; 
                    [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
                    glm_cellinfo.WN_STA = STAandROI.STA;
                    glm_cellinfo.slave_centercoord = center_coord;
                    clear master_idx stafit_centercoord slvdim sd center_coord sd

                    %% Execute GLM FITTING
                    
                    display(sprintf('### running: %s %s %s ###', GLMType.fit_type, expname, cell_savename))
                    display(sprintf('### fit type: %s ###', GLMType.fitname))
                    % EXECUTE THE CORRECT FITTING
                    if isfield(GLMType, 'DoubleOpt') && GLMType.DoubleOpt
                        if isfield(GLMType,'DoubleOpt_Manual') && GLMType.DoubleOpt_Manual
                            [fittedGLM, manual_search] = glm_execute_DoubleOpt_Manual(GLMType, spikesconcat, concat_fitmovie, inputstats, glm_cellinfo);
                        else
                            [fittedGLM] = glm_execute_DoubleOpt(GLMType, spikesconcat, concat_fitmovie, inputstats, glm_cellinfo);
                        end
                    else
                        [fittedGLM]     = glm_execute(GLMType, spikesconcat, concat_fitmovie, inputstats, glm_cellinfo);            
                    end
                    
                    if strcmp(GLMType.stimfilter_mode, 'rk2-ConductanceBased') || strcmp(GLMType.stimfilter_mode, 'fixedSP-ConductanceBased')
                        xvalperformance = hack_CB_eval_xvalperformance(fittedGLM, StimulusPars.slv, organizedspikes,testmovie,inputstats);
                        fittedGLM.xvalperformance  = xvalperformance; 
                        fittedGLM.d_save           = Dirs.fittedGLM_savedir;
                        eval(sprintf('save %s/%s.mat fittedGLM', Dirs.fittedGLM_savedir, glm_cellinfo.cell_savename));
                        %
                        printname = sprintf('%s/DiagPlots_%s', Dirs.fittedGLM_savedir,fittedGLM.cellinfo.cell_savename);
                        hack_CB_printglmfit(fittedGLM,printname)
                    else
                        % EVALUATE THE FITS: SAVE AND PLOT
                        xvalperformance = eval_xvalperformance_NEW(fittedGLM, StimulusPars.slv, organizedspikes,testmovie,inputstats);
                        fittedGLM.xvalperformance  = xvalperformance; 
                        fittedGLM.d_save           = Dirs.fittedGLM_savedir;
                        eval(sprintf('save %s/%s.mat fittedGLM', Dirs.fittedGLM_savedir, glm_cellinfo.cell_savename));
                        %
                        printname = sprintf('%s/DiagPlots_%s', Dirs.fittedGLM_savedir,fittedGLM.cellinfo.cell_savename);
                        printglmfit(fittedGLM,printname)
                    end
                    
                    % EXTRA PARAM CHECK FOR DOUBLE OPT MANUAL 
                    %{
                    if isfield(GLMType,'DoubleOpt_Manual') && GLMType.DoubleOpt_Manual
                        paramdir = sprintf('%s/ParamLandscape', Dirs.fittedGLM_savedir);
                        if ~exist(paramdir, 'dir'), mkdir(paramdir); end
                        eval(sprintf('save %s/%s.mat manual_search', paramdir, glm_cellinfo.cell_savename));
                        param_printname = sprintf('%s/ParamLandscape_%s', paramdir,fittedGLM.cellinfo.cell_savename);
                        printparamvar(fittedGLM,manual_search,param_printname)
                    end
                    %}
                end
            end
        end
    end
end

end