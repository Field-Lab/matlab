clear

%% Classification run
Analysis_Path = '/Volumes/Analysis/2015-05-27-3/data001-data005';
datarun_class = load_data([Analysis_Path '/data001/data001'], struct('load_neurons', 0, 'load_params', 1));
model_fit_data = load_data([Analysis_Path '/data002/data002'], struct('load_neurons', 1, 'load_params', 0));
dsave = '/Volumes/Lab/Users/Nora/GLMFits/2015-05-27-3';

%% ONLY CHANGE THINGS HERE
cell_spec = get_cell_ids(datarun_class,'On Parasol'); % cell ids to fit
%cell_spec = [202];
convergence = 1; % fraction of data to use

%% Don't change these
if convergence < 1; dsave = [dsave '_Conv_' num2str(convergence)]; end
mkdir(dsave)

block_starts = [1 503 1485 1987 2969 3471 4453 4955 5937 6439];
blocks{1} = block_starts(1:2:end);
blocks{2} = block_starts(2:2:end);
monitor_refresh = 120;
visual_check = 1;


%%
% block frames is a vector of the different block lengths,
% ex = [1200*3 1200] aka even fitting blocks first and odd testing blocks
% second. Weird order is so that frames(block i) = block_frames(mod(i,2)+1)
%
% Output:
% preppeddata structure, including start_time which has the starts of the
% blocks and p which are the optional parameters you included
% It can also include: fitspikes, testspikes, fitmovie, testmovie, etc
% which depend on which options you choose. They will be in cells and each
% cell corresponds to a block of data.
%
% cell_spec
% Setting cell_specification will return the organized spikes for those
% cells. Details for what cell_specification should be are detailed in the
% get_cell_indices function

%% Find Block Times
clear prepped_data
monitor_refresh = 120;
visual_check = 1;

for i_stim = 2
    trial_idx = 1:10;
    tic;
    for i_block_trigger = blocks{i_stim}
        
        % load up the triggers and find triggers per block
        try
            triggers = model_fit_data.triggers(i_block_trigger+(1:1000)-1);
        catch
            triggers = model_fit_data.triggers(i_block_trigger:end);
        end
        block_frames = i_stim*[3600 1200];
        n_repeats = 10;
        
        % calculate some things
        repeat_only = false;
        n_blocks = n_repeats*2;
        testmovie_only = 0;
        triggers_per_block = block_frames/100;
        
        % Initialize and define things
        start_trigger = zeros(n_blocks,1);
        start_trigger(1) = 1;
        
        % Find the trigger numbers of the block starts
        for i = 1:n_blocks
            end_trigger = start_trigger(i) +  triggers_per_block(mod(i,2)+1);
            start_trigger(i+1) = end_trigger + 1;
        end
        start_trigger = start_trigger(start_trigger <= length(triggers));
        n_blocks = length(start_trigger);
        
        % visual check
        if 0
            figure;
            plot(diff(triggers))
            hold on
            plot(start_trigger, 1/1.2*ones(length(start_trigger)), '*')
            title('Identifying the Block Starts')
        end
        
        % find the start times
        start_time = triggers(start_trigger);
        prepped_data.start_time = start_time;
        
        % Find and organize the spikes
        % Organize the cell spikes
        if cell_spec
            if isstruct(datarun_class)
                cids = get_cell_indices(datarun_class, cell_spec);
            else
                cids = get_cell_indices(model_fit_data, cell_spec);
            end
            spikes = cell(n_blocks,length(cids)); % initialize spike cell array
            for i_cell = 1:length(cids)
                spikes_temp = model_fit_data.spikes{cids(i_cell)};
                for i_trigger = 1:n_blocks
                    block_length = block_frames(mod(i_trigger,2)+1)/monitor_refresh;
                    spikes{i_trigger, i_cell} = spikes_temp( (spikes_temp > start_time(i_trigger)) & (spikes_temp < (start_time(i_trigger) + block_length)) ) - start_time(i_trigger);
                end
            end
            if repeat_only
                prepped_data.testspikes = spikes;
            else
                prepped_data.fitspikes(trial_idx,:) = spikes(2:2:20, :);
                prepped_data.testspikes(trial_idx,:) = spikes(1:2:19, :);
                trial_idx = trial_idx+10;
            end
            clear spikes
        end
    end
    
    if 0
        figure;
        hold on
        for i = 1:size(prepped_data.testspikes, 1)
            plot(prepped_data.testspikes{i, 1}, i*ones(length(prepped_data.testspikes{i, 1})), 'k.')
        end
        title('Checking Cell Repeats')
    end
    
    disp(['Trigger and spike finding took ' num2str(toc) ' seconds.'])
    
    % Load up and organize the stimulus
    tic;
    if i_stim == 1
        WN_stim = interleaved_data_prep(model_fit_data, [3600 1200], 50, 'stimulus_name', 'BW-8-1');
        prepped_data.testmovie = WN_stim.testmovie;
        save_name = '';
    else
        load('/Volumes/Lab/Users/Nora/downsampledNSbrownian.mat');
        fitmovie = fitmovie(:,:,2401:end);
        stimlength = round(size(fitmovie,3)*convergence);
        if convergence < 1
            fitmovie = fitmovie(:,:,1:stimlength);
        end
            
        test = load('/Volumes/Lab/Users/Nora/downsampledNSbrownian_testA.mat');
        prepped_data.testmovie = test.fitmovie;
        save_name = 'NSEM';
    end   
    
    % Concat WN movie
    if i_stim == 1;
        grey_buffer = 0;
        n_blocks = size(WN_stim.fitmovie,1);
        block_size = size(WN_stim.fitmovie{1});
        block_length = block_size(3);
        grey_frames = 60;
        if grey_buffer
            block_length = block_length+grey_frames;
            movie_mean = mean(WN_stim.fitmovie{1}(:));
            grey_movie = movie_mean*ones(block_size(1), block_size(2), grey_frames);
        end % one STA length
        block_size(3) = block_length*n_blocks; %movie size
        fitmovie = zeros(block_size,'uint8');
        
        % concatenate the movie
        idx = 1:block_length; 
        for i_block = 1:n_blocks
            if ~grey_buffer
                fitmovie(:,:,idx) = WN_stim.fitmovie{i_block};
            else
                fitmovie(:,:,idx) = cat(3, WN_stim.fitmovie{i_block}, grey_movie);
            end
            idx = idx+block_length;
        end
        if convergence < 1; fitmovie = fitmovie(:,:,1:(block_length*n_blocks*convergence)); end
    end
    disp(['Stimulus organization took ' num2str(toc) ' seconds.']); tic
      
    % fit each cell
    n_cells = size(prepped_data.fitspikes,2);
    block_length = i_stim*[3600];
    n_blocks = size(prepped_data.fitspikes,1);
    for i_cell=11:n_cells
        cell_savename = num2str(cell_spec(i_cell));
        fitspikes = [];
        for i_block = 1:n_blocks
            t_block_start = block_length*(i_block - 1)/monitor_refresh;
            fitspikes = [fitspikes; prepped_data.fitspikes{i_block,i_cell}+t_block_start];
        end
        fitspikes = fitspikes(fitspikes < stimlength/monitor_refresh);
        
        if i_stim == 1
            close all
            [STA, ~] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
            center = round([40 - datarun_class.vision.sta_fits{cids(i_cell)}.mean(2) datarun_class.vision.sta_fits{cids(i_cell)}.mean(1)]);
            hold on; plot(center(1), center(2), 'r*')
        else
             eval(sprintf('load %s/%s.mat fittedGLM', dsave, cell_savename));
             STA = fittedGLM.STA;
             center = fittedGLM.center;
             clear fittedGLM
        end
    
        fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA', STA, 'monitor_refresh', monitor_refresh);
        fittedGLM.xvalperformance = glm_predict(fittedGLM, prepped_data.testmovie, 'testspikes', prepped_data.testspikes(:,i_cell));
        fittedGLM.STA = STA;
        fittedGLM.center = center;
        
        save([dsave '/' cell_savename save_name '.mat'],'fittedGLM', '-v7.3');
        close all
        plotfilters(fittedGLM);
        set(gcf, 'Position', [100 100 800 250])
        exportfig(gcf, [dsave '/' cell_savename '_' save_name 'filters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
        
        plotrasters(fittedGLM.xvalperformance, fittedGLM);
        exportfig(gcf, [dsave '/' cell_savename '_' save_name 'rasters'], 'Bounds', 'loose', 'Color', 'rgb');
    end
end
