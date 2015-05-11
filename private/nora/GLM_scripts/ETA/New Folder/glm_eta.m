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
troubleshoot.name    = 'singleopt';

%  LOOP THROUGH DATA SETS

BD = NSEM_BaseDirectories;

exptests = [1];
cellselectiontype = 'shortlist';
troubleshoot.plotdir = BD.GLM_troubleshootplots;

%%

count = 0;
for i_exp = exptests
   
    expnumber = i_exp;
    [exp_nm,cells,~]  = cell_list( expnumber, cellselectiontype);
    [StimulusPars, ~, ~, datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    inputs.exp_nm    = exp_nm;
    inputs.map_type  = GLMType.map_type;
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs;
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));
    display(sprintf('Save Directory :  %s', d_save));
    
    % Load movie
    clear Main_SolPars Other_SolParss
    [testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
    
    % Loop through cells
    for i_cell = 1:length(cells)
        clear glm_cellstruct
        cid = cells{i_cell};
        [~ , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);

        % Load the fit GLM
        try
            load([d_save '/' cell_savename '.mat'])
            disp('Fitted GLM Loaded');
            
            % Calculate the PSTH
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
            
            % maybe saccade input has to do with change in stimulus in nearby zone?
            if 1
                sta = fittedGLM.linearfilters.Stimulus.space_rk1;
                x_coord = fittedGLM.linearfilters.Stimulus.x_coord;
                y_coord = fittedGLM.linearfilters.Stimulus.y_coord;
                glm_error = PSTH(1,:) - PSTH(2,:);
                buffer = 1000;
                err_count = 0;
                for i = 1:28
                    bin_number = round(i / dt);
                    frame_number = round(i * 120);
                    idx = (bin_number-buffer):(bin_number+buffer);
                    temp = max(glm_error(idx));
                    temp_min = min(glm_error(idx));
                    if abs(temp_min) > temp
                        temp = temp_min;
                    end
                    if abs(temp) > 200
                        err_count = err_count + 1;
                        big_errors{i_cell, err_count}.saccade = i;
                        big_errors{i_cell, err_count}.before = testmovie{1}.matrix(x_coord, y_coord, frame_number-50);
                        big_errors{i_cell, err_count}.after = testmovie{1}.matrix(x_coord, y_coord, frame_number+50);
                    end
                    SI(i_cell, i) = temp;
                end
                figure;
                for i=1:err_count
                    subplot(err_count, 2, 2*i-1)
                    image(big_errors{i_cell, i}.before)
                    axis image
                    title(['Saccade' num2str(big_errors{i_cell, i}.saccade)])
                    subplot(err_count, 2, 2*i)
                    image(big_errors{i_cell, i}.after)
                    axis image
                end
            end
            
        catch
        end
    end
    
end


%%

for i_cell=1:9
    err_count = 0;
    figure;
    for i_err = 1:10
        if isstruct(big_errors{i_cell, i_err})
            err_count = err_count+1;
        end
    end
    for i_err=1:err_count
        subplot(1,err_count, i_err)
        image(big_errors{i_cell, i_err}.before - big_errors{i_cell, i_err}.after)
        axis image
    end
end


