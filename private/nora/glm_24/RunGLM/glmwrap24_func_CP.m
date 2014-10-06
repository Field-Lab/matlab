% The inner loop should be a straight copy of glmwrap24
% Develop in glmwrap24,  transfer over here for other
% AKHeitman 2014-03-27
% AKHeitman 2014-04-07  -emphasize role as parameter independet loading
%                       -setting directories
%                       -dicatating parameters
%                       -no computations
% AKHeitman 2014-06-08  -make functional version of glmwrap
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



%%  LOOP THROUGH DATA SETS


function glmwrap24_func_CP(GLMType,BaseDirectories, exptests,cellselectiontype)
BD = BaseDirectories;
GLMType.func_sname = 'glmwrap24_func';
GLMType.fullmfilename =mfilename('fullpath'); 

for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    cells
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    
        % NBCoupling 06-12-14
    %lets start using multiple pairings, still manual
    % cells={5086, 841,1471}; % cells to fit
    % cell_pairs{1}=[1471 5161]; % cells to pair with the first cell
    % cell_pairs{2}=[1426]; % cells to pair with the second cell
    % cell_pairs{3}=[5086, 1276]; % etc
    if GLMType.CouplingFilters==true
        CTYPE = {'On-Parasol','Off-Parasol'}; %
        for n_ctype = 1:length(CTYPE)
            [~, parasol_cp.ids{n_ctype}] = celltype_id_AH(CTYPE{n_ctype}, datarun_slv.cell_types);
        end
        parasol_cp.indices{1}=get_cell_indices(datarun_mas,parasol_cp.ids{1});
        parasol_cp.indices{2}=get_cell_indices(datarun_mas,parasol_cp.ids{2});
        parasol_cp.NumCells{1}=length(parasol_cp.ids{1});
        parasol_cp.NumCells{2}=length(parasol_cp.ids{2});
    end
    % end NBCoupling
    
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
    [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
    [testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
    GLMType.fitmoviefile = origmatfile;
    clear origmatfile
    %uint8 form to keep data low


    concat_fitmovie      = concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv);
    fitmoviestats.mean   = mean(concat_fitmovie(:));
    fitmoviestats.minval =  min(concat_fitmovie(:));
    fitmoviestats.maxval =  max(concat_fitmovie(:));
    clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
    %}
 
    %% Load Cell Specific Elements   Spikes and STA
    inputs.exp_nm       = exp_nm; 
    inputs.map_type     = GLMType.map_type; 
    DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs); 
    inputs.stim_type    = GLMType.fit_type;
    DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs); 
    clear inputs

    for i_cell = 1:length(cells)

        clear glm_cellstruct
        cid = cells{i_cell};
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);  
        if ~exist(sprintf('%s/%s.mat', d_save,cell_savename))
            glm_cellinfo.cid           = cid;
            glm_cellinfo.exp_nm        = exp_nm;
            glm_cellinfo.celltype      = celltype;
            glm_cellinfo.cell_savename = cell_savename;
            glm_cellinfo.fitname       = GLMType.fitname
            glm_cellinfo.computedtstim = StimulusPars.slv.computedtstim;
            %%% only really matters when binning spikes for fitting %%
            %%% small differences from 1/120 will start adding up over the long run
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Load the center_coord  (do Param dependent ROI manipulation later
            master_idx         = find(datarun_mas.cell_ids == cid);
            stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
            stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
            slvdim.height      = StimulusPars.slv.height; slvdim.width = StimulusPars.slv.width; 
            [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
            glm_cellinfo.slave_centercoord = center_coord;
            clear master_idx stafit_centercoord slvdim sd
            
            
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
            spikesconcat.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
            eval(sprintf('load %s/STAandROI_%s.mat STAandROI', DirPars.WN_STAdir, cell_savename));
            glm_cellinfo.WN_STA = STAandROI.STA;
            clear cell_savename
            
            % NBCoupling 05-28-14
            % Make a neighbor spikes loop here! %
            % Eventually figure out a way to not have to load spikes twice
            if GLMType.CouplingFilters
                n_couplings=length(glm_cellinfo.pairs); % number of cells to couple to
                % loading the neighboring spikes to neighborspikes.home
                for j=1:n_couplings
                    [~ , neighborspikes.pair_savename{j}, ~]  = findcelltype(glm_cellinfo.pairs(j), datarun_mas.cell_types);
                    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, neighborspikes.pair_savename{j}));
                    neighborspikes.home{j} = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                end
            else
                % if there's no coupling, just set this to zero
                neighborspikes=0;
            end
            % end NBCoupling

            % Prepare the stim
            % DO SOMETHING ABOUT COMPUTED TSTIM!!! 
            
            %% Execute the correct GLM
            tic
            if isfield(GLMType, 'DoubleOpt') && GLMType.DoubleOpt
                [fittedGLM] =    glm_execute_DoubleOpt(GLMType, spikesconcat, concat_fitmovie, glm_cellinfo);
            else
                [fittedGLM]     = glm_execute_CP(GLMType, spikesconcat, concat_fitmovie, glm_cellinfo);            
            end
            toc
            
            xvalperformance = eval_xvalperformance_NEW(fittedGLM, StimulusPars.slv, organizedspikes,testmovie);
            fittedGLM.xvalperformance  = xvalperformance; 
            fittedGLM.d_save           = d_save;
            eval(sprintf('save %s/%s.mat fittedGLM', d_save, glm_cellinfo.cell_savename));
            printname = sprintf('%s/DiagPlots_%s', d_save,fittedGLM.cellinfo.cell_savename);
            printglmfit(fittedGLM,printname)
       % NB 06-11-2014
        else 
            error('Previous results still in directory')
        end
        
         
    end
    
end
