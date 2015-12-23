% combine chunks into a single datarun, save as mat file


tic

% parameters
save_what = 'sta';
%save_what = 'rf';

frames_to_save = 5;


% choose dataset folder
clear sources
switch tsc
    case 1 % blueberry, 20 min chunks
        master = '2008-08-26-2/data001/data001';
        root_dir = '2008-08-26-2/data001/pieces/';
        for ss=1:8; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        cell_spec = 'all';
        keep_chunks = [1 2 3];
    case 3 % kiwi, 15 min chunks
        master = '2008-05-13-3/data006/data006';
        root_dir = '2008-05-13-3/data006/pieces/';
        for ss=1:9; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        cell_spec = 'all';
        keep_chunks = [1 2 3];
        
    case 7 % cherimoya, 10 min chunks
        master = '2008-03-25-3/data002/data002';
        root_dir = '2008-03-25-3/data002/pieces/';
        for ss=1:7; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        cell_spec = 'all';
        keep_chunks = [5:8];
end


% initialize with load master dataset info
datarun_sum = load_data(master);
datarun_sum = load_sta(datarun_sum,'load_sta',[],'verbose',1);
datarun_sum = load_params(datarun_sum,struct('verbose',1));

% load basic info from each datarun
datarun = load_data({sources{keep_chunks}});
for dd=1:length(datarun)
    datarun{dd} = load_sta(datarun{dd},'load_sta',[],'verbose',1);
    datarun{dd} = load_neurons(datarun{dd});
end

% sum up STAs from each datarun, store in datarun_sum

% identify number of frames
num_frames = size(get_sta(datarun{1},datarun{1}.cell_ids(1)),4);

% go throuch each cell
cell_indices = get_cell_indices(datarun_sum,cell_spec);
for cell_index=cell_indices
    fprintf('.')
    cell_id = datarun_sum.cell_ids(cell_index);

    % get data from each chunk
    
    % initialize variables to store data
    sta_temp = zeros(datarun_sum.stimulus.field_height,datarun_sum.stimulus.field_width,3,num_frames);
    spikes_temp = [];
    
    % go through each datarun, note its part of the data
    for dd =1:length(datarun)
        % check if the cell id is found in this datarun
        if any(datarun{dd}.cell_ids==cell_id)
            
            ci = get_cell_indices(datarun{dd},cell_id);
            
            % get spikes
            spikes_temp = [spikes_temp; datarun{dd}.spikes{ci}];
            
            % get STA
            sta_partial = get_sta(datarun{dd},cell_id);
            % multiply by number of spikes
            sta_partial = sta_partial*length(datarun{dd}.spikes{ci});
            % add it in
            sta_temp = sta_temp + sta_partial;
            
        end
    end
    
    
    % save to datarun_sum
    
    % save spikes
    datarun_sum.spikes{cell_index} = sort(spikes_temp);
    
    % save STA or RF
    switch save_what
        case 'sta' % save STA
            datarun_sum.stas.stas{cell_index} = sta_temp(:,:,:,end-frames_to_save+1:end) / ...
                length(datarun_sum.spikes{cell_index});
        case 'rf' % save RF
            % convert to RF
            % get sig_stixels
            ss_params = struct;
            ss_params.strength = {'inner or',[0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};
            sig_stixels = significant_stixels(sta_temp,ss_params);
            % get RF
            rf = rf_from_sta(sta_temp,'sig_stixels',sig_stixels);
            % store in datarun_sum
            datarun_sum.stas.rfs{cell_index} = rf;
        otherwise
            error('must save SOMETHING!')
    end
    
end
fprintf('\n')
toc

% compute typical STA summaries
datarun = get_sta_summaries(datarun,cell_spec,'keep_rfs',0,'verbose',1,'fig_or_axes',3);


