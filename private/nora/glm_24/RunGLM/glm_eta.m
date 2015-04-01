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
clear; close all;  clc

%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE!

% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS

% SETUP cells and experiments, the TYPE of GLM (GLMType)

GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean';
GLMType.nullpoint = 'mean';
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false;
GLMType.specialchange = false;
GLMType.specialchange_name = 'extra_coupling';
GLMType.CBP=false;

GLMType.stimfilter_mode = 'rk1';
%GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = false;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';

GLMType.CONVEX = true;
GLMType.DoubleOpt = false;
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}

GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = true;
GLMType.Subunits = false;
GLMType.Saccades = false;
GLMType.color = false;

% GLMType.fixed_spatialfilter = true;
% NBCoupling 06-12-2014
GLMType.func_sname = 'glmwrap24_CP';
GLMType.fullmfilename =mfilename('fullpath');
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);
troubleshoot.doit    = false;
%troubleshoot.plotdir = '/Users/akheitman/Matlab_code/troubleshooting_plots'
% NBCoupling
%troubleshoot.plotdir=BD.GLM_troubleshootplots;
troubleshoot.name    = 'singleopt';

%  LOOP THROUGH DATA SETS

BD = NSEM_BaseDirectories;

exptests = [1 2 3 4 ];
cellselectiontype = 'shortlist';
troubleshoot.plotdir = BD.GLM_troubleshootplots
%%
count = 0;
for i_exp = exptests
    %%
    expnumber = i_exp;
    [exp_nm,cells,~]  = cell_list( expnumber, cellselectiontype);
    [StimulusPars, ~, ~, datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    
    % NBCoupling 06-12-14
    %     if GLMType.CouplingFilters==true
    %         CTYPE = {'On-Parasol','Off-Parasol'};
    %         for n_ctype = 1:length(CTYPE)
    %             [~, parasol_cp.ids{n_ctype}] = celltype_id_AH(CTYPE{n_ctype}, datarun_slv.cell_types);
    %         end
    %         parasol_cp.indices{1}=get_cell_indices(datarun_mas,parasol_cp.ids{1});
    %         parasol_cp.indices{2}=get_cell_indices(datarun_mas,parasol_cp.ids{2});
    %         parasol_cp.NumCells{1}=length(parasol_cp.ids{1});
    %         parasol_cp.NumCells{2}=length(parasol_cp.ids{2});
    %     end
    % end NBCoupling
    
    %     %%%%  Shorten Block count if using Debug
    %     if GLMType.debug
    %         StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:10);
    %     end
    %     clear boolean_debug map_type fit_type shead_cellID expname
    %
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
    
    inputs.exp_nm    = exp_nm;
    inputs.map_type  = GLMType.map_type;
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
    
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs;
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));
    display(sprintf('Save Directory :  %s', d_save));
    %     if ~exist(d_save), mkdir(d_save); end
    GLMType.d_save = d_save;
    
    %% Load Movie and Concatenate the Fitting Section
    clear Main_SolPars Other_SolParss
    %%% Load Stimulus   -- insert more frame cutting here!
    %     [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
    [testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
    %     GLMType.fitmoviefile = origmatfile;
    %     clear origmatfile
    %uint8 form to keep data low
    
    
    %     concat_fitmovie      = concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv);
    %     fitmoviestats.mean   = mean(concat_fitmovie(:));
    %     fitmoviestats.minval =  min(concat_fitmovie(:));
    %     fitmoviestats.maxval =  max(concat_fitmovie(:));
    %     clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
    %     %}
    
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
        %             % NBCoupling 06-10-2014
        %             if GLMType.CouplingFilters==true
        %                 glm_cellinfo.pairs=pick_neighbor_cells(stafit_centercoord, parasol_cp, datarun_mas.vision.sta_fits);
        %             else
        %                 glm_cellinfo.pairs=0;
        %             end
        %             % end NBCoupling
        %             clear master_idx stafit_centercoord slvdim sd
        %
        %             eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
        %             spikesconcat.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
        %             eval(sprintf('load %s/STAandROI_%s.mat STAandROI', DirPars.WN_STAdir, cell_savename));
        %             glm_cellinfo.WN_STA = STAandROI.STA;
        %             clear cell_savename
        
        %
        %             % NBCoupling 05-28-14
        %             % Make a neighbor spikes loop here! %
        %             % Eventually figure out a way to not have to load spikes twice
        %             cell_organizedspikes=organizedspikes;
        %             if GLMType.CouplingFilters
        %                 n_couplings=length(glm_cellinfo.pairs); % number of cells to couple to
        %                 % loading the neighboring spikes to neighborspikes.home
        %                 for j=1:n_couplings
        %                     [~ , glm_cellinfo.pair_savename{j}, ~]  = findcelltype(glm_cellinfo.pairs(j), datarun_mas.cell_types);
        %                     eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, glm_cellinfo.pair_savename{j}));
        %                     neighborspikes.home{j} = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
        %                     neighbor_organizedspikes{j}=organizedspikes;
        %                 end
        %                 if GLMType.Saccades
        %                    neighborspikes.home{n_couplings+1}=(0:1:(size(concat_fitmovie,3)/120))';
        %                    neighbor_organizedspikes{n_couplings+1}.block.t_sp_withinblock=repmat({0:StimulusPars.slv.seconds_perstaticblock;0:StimulusPars.slv.seconds_perstaticblock},StimulusPars.slv.n_rep,1);
        %                 end
        %             else
        %                 if GLMType.Saccades
        %                     neighborspikes.home{1}=(0:1:(size(concat_fitmovie,3)/120))';
        %                     neighbor_organizedspikes{1}.block.t_sp_withinblock=repmat({0:StimulusPars.slv.seconds_perstaticblock;0:StimulusPars.slv.seconds_perstaticblock},StimulusPars.slv.n_rep,1);
        %                 else
        %                     neighborspikes=0;neighbor_organizedspikes=0;
        %                 end
        %             end
        %             % end NBCoupling
        %
        %             % Prepare the stim
        %             % DO SOMETHING ABOUT COMPUTED TSTIM!!!
        %% Load the correct GLM
        
        try
            load([d_save '/' glm_cellinfo.cell_savename '.mat'])
            disp('Fitted GLM Loaded');
            
            %% Plot the PSTH
            convolve=100;
            rec_rast = fittedGLM.xvalperformance.rasters.recorded;
            sim_rast = fittedGLM.xvalperformance.rasters.glm_sim;
            trials = size(rec_rast,1);
            for i=1:trials
                PSTH_rec(i,:)=conv(rec_rast(i,:),ones(convolve,1),'same');
                PSTH_sim(i,:)=conv(sim_rast(i,:),ones(convolve,1),'same');
            end
            PSTH = zeros(2, size(PSTH_rec,2));
            PSTH(1,:) = sum(PSTH_rec);
            PSTH(2,:) = sum(PSTH_sim);
            dt = fittedGLM.t_bin;
            time = dt:dt:dt*size(PSTH,2);
            disp('PSTH calculated');
            
            % Is the amount of error related to the strength of the signal?
            plot(PSTH(1,:), abs(PSTH(1,:) - PSTH(2,:)),'.');
            pause()
           
            
            %% maybe saccade input has to do with change in stimulus in nearby zone?
            
            if 0
                sta = fittedGLM.linearfilters.Stimulus.space_rk1;
                x_coord = fittedGLM.linearfilters.Stimulus.x_coord;
                y_coord = fittedGLM.linearfilters.Stimulus.y_coord;
                glm_error = PSTH(1,:) - PSTH(2,:);
                buffer = 1000;
                figure;
                hold on
                for i = 1:28
                    bin_number = round(i / dt);
                    frame_number = round(i * 120);
                    before_stim = sum(sum(double(testmovie{1}.matrix(x_coord, y_coord, frame_number-100))*sta));
                    after_stim = sum(sum(double(testmovie{1}.matrix(x_coord, y_coord, frame_number+100))*sta));
                    idx = (bin_number-buffer):(bin_number+buffer);
                    plot(glm_error(idx));
                    SI(1,i) = max(abs(glm_error(idx)));
                    SI(2,i) = abs(after_stim - before_stim);
                end
                plot([buffer+1 buffer+1], [-500, 500], 'k')
                xlabel('Bin')
                ylabel('GLM PSTH - True PSTH (arb. units)')
                xlim([0 2*buffer+1])
                text(1010, 300, 'Time of Saccade')
                title(cell_savename)
                hold off
                pause()
            end
            
                        
            %% are most of the errors even by the saccades?
            
            if 0
                
                count = count + 1;
                sta = fittedGLM.linearfilters.Stimulus.space_rk1;
                x_coord = fittedGLM.linearfilters.Stimulus.x_coord;
                y_coord = fittedGLM.linearfilters.Stimulus.y_coord;
                glm_error = PSTH(1,:) - PSTH(2,:);
                buffer = 100;
                saccade_error = 0;
                n_bins = 0;
                for i = 1:28
                    bin_number = round(i / dt);
                    idx = (bin_number-1):(bin_number+buffer);
                    saccade_error = saccade_error + sum(abs(glm_error(idx)));
                    n_bins = n_bins + length(idx);
                end
                s_e(1,count) = saccade_error/n_bins;
                s_e(2,count) = sum(abs(glm_error))/length(glm_error);
                disp(saccade_error/sum(abs(glm_error)))
                
            end

            %% Calculate the ETA
            
            if 0
            tic;
            glm_error = PSTH(1,:) - PSTH(2,:);
            eta = zeros(80,40,30);
            control = zeros(80,40,30);
            for i=1:length(glm_error)
                frame = ceil(i * dt * 120)+120;
                if abs(glm_error(i))>100
                    temp = double(testmovie{1}.matrix(:,:,(frame-29):frame));
                    eta = eta + temp;
                end
                temp = double(testmovie{1}.matrix(:,:,(frame-29):frame));
                control = control + temp;
            end
            disp('ETA calculated');
            toc
            
            movie_file = VideoWriter(['/Users/Nora/Desktop/ETAs/' cell_savename '.avi']);
            open(movie_file)
            for i =1:30
                subplot(2,1,1)
                imagesc(eta(:,:,i)');
                hold on; plot(stafit_centercoord(1),stafit_centercoord(2),'o y'); hold off
                colormap gray
                axis image
                axis off
                title('Error Triggered Average')
                subplot(2,1,2)
                imagesc(control(:,:,i)');
                hold on; plot(stafit_centercoord(1),stafit_centercoord(2),'oy'); hold off
                colormap gray
                axis image
                axis off
                title('Control Average')
                pause()
                writeVideo(movie_file, getframe(gcf))
                writeVideo(movie_file, getframe(gcf))
            end
            close(movie_file)
            end
            
        catch
        end
    end
    
end


%%

