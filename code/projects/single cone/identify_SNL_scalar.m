% for each cell, identify scale factors to fit static nonlinearity
% fortunately, each the scale factors depend very little on whether extra noise stixels are added to the RF

% assume datarun is already loaded


% PARAMETERS

% which cells
cell_spec = {1,2,3,4,5,6}; cell_spec = get_cell_ids(datarun,cell_spec); cell_spec = cell_spec(1:10);
%cell_spec = 'all';
%cell_spec = cell_ids;

% first stimulus to use
first_stim = 1;

% how many stimuli to use
num_stims = 500;

% frames to use
num_frames = 3;



% INITIALIZE STUFF

% get stimulus movie and refresh time
if ~exist('refresh_time','var')
    movie = edu.ucsc.neurobiology.vision.matlab.Matlab.computeMovie(datarun.stimulus.xml_file, datarun.triggers(1:100)*20000);
    refresh_time = movie.getRefreshTime/1000;
end

% note stimulus size
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% identify which stim frames to use
stim_start = first_stim + num_frames - 1; % the first few can't be used, because the STRF occupies several frames
stim_end =  first_stim + num_stims - 1 ;




% GET STRFs AND SPIKE TIMES
if 1

    % initialize storage variables
    all_spike_times = sparse(zeros(num_stims-num_frames+1,length(cell_indices)));
    all_sig_stixels = sparse(false(field_width*field_height,length(cell_indices)));
    all_strfs = cell(length(cell_indices),1);

    % cycle through each cell
    for cc = 1:length(cell_indices)
        fprintf('.')

        cell_index = cell_indices(cc);
        cell_id = datarun.cell_ids(cell_index);


        % get spike times

        % compute spike rate at all times
        spike_rate_ = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.duration);
        % store spikes in the relevant region
        all_spike_times(:,cc) = spike_rate_(stim_start:stim_end);


        % get sig stixels
        switch 1
            case 1 % marks
                sig_stixels = datarun.stas.marks{cell_index};
            case 2 % identify sig stixels now
                sig_stixels = significant_stixels(get_sta(datarun,cell_id));
        end
        % reshape
        sig_stixels = reshape(sig_stixels,[],1);
        % store
        all_sig_stixels(:,cc) = sig_stixels;


        % get STRF

        % get STA
        sta = get_sta(datarun,cell_id);
        % pare to relevant frames
        sta = sta(:,:,:,end-num_frames+1:end);
        % reshape, so that first dim is space, second dim is time-color
        sta_r = reshape(sta,[],size(sta,4));
        % pare to relevant region
        strf = reshape(sta_r(repmat(full(all_sig_stixels(:,cc)),3,1),:),[],1);
        % store
        all_strfs{cc} = strf;

    end
    fprintf('\n')

end







% GET GENERATOR SIGNAL VALUES
if 1

    % initialize storage variables
    all_gen_signal = zeros(num_stims-num_frames+1,length(cell_indices));
    all_fit_params = zeros(2,length(cell_indices));
    

    % load up initial stimuli
    % initialize stimulus
    clear stims
    for ss=first_stim:first_stim+num_frames-2
        % get new frame
        STAFrame = movie.getFrame(ss-1);
        new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

        % assign this frame to the current position
        stims(:,ss+2-first_stim) = reshape(new_frame,[],1);
    end

    % cycle through each stimulus
    for ss = stim_start:stim_end
        
        % LOAD STIMULUS
        
        % get new frame
        STAFrame = movie.getFrame(ss-1);
        new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

        % shift current frames over
        stims(:,1:end-1) = stims(:,2:end);
        % assign new frame to the final position
        stims(:,end) = reshape(new_frame,[],1);
        
        %disp(stims(1:6,:));pause
        
        
        % COMPUTE GENERATOR SIGNAL VALUE FOR EACH CELL

       
        % get generator signal values for each cell
        for cc = 1:length(cell_indices)
            % get STRF
            strf = all_strfs{cc};

            % get stimulus in relevant region
            stim = reshape(stims(repmat(full(all_sig_stixels(:,cc)),3,1),:),[],1);

            % take dot product
            stim_num = ss - first_stim - num_frames + 2;
            all_gen_signal(stim_num,cc) = stim'*strf;
        end
        
        
    end
    
    
    

    % FIT PARAMETERS AND STORE THEM
    
    % cycle through each cell
    for cc = 1:length(cell_indices)
        
        % get spike times
        %spike_times = find(all_spike_times(:,cc));
        

        % get spike times
        % if there are N spikes per bin, list that bin N times
        spike_times = [];
        for nn = 1:full(max(all_spike_times(:,cc)))
            spike_times = [spike_times; find( all_spike_times(:,cc)>(nn-1) )];
        end

        % define function that gives likelihood
        L = @(x)qd_identify_SNL_scalars_fn(x,spike_times,all_gen_signal(:,cc));

        % turn gradient on
        options = optimset('GradObj','on');

        % identify best a and b
        all_fit_params(:,cc) = fminunc(L,[0 0],options)';
        
    end
    
end




% PLOT RESULTS
if 1
    figure(2);clf;hold on
    cols = 'rgbk';
    ct = find_cell_types(datarun,datarun.cell_ids(cell_indices));
    for tt=1:4
        plot(all_fit_params(1,ct==tt),exp(all_fit_params(2,ct==tt)),'.','Color',cols(tt))
    end
end
