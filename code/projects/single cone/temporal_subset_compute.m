% Identify which temporal subset of a single cone RF run has the least motion
%
% STAs should be computed for 10-20 minute chunks of the run, and stored in this format:
%
% 2008-08-26-2/data001/pieces/part1/part1.sta
% 2008-08-26-2/data001/pieces/part2/part2.sta
%  ...
%
% This script then combines adjacent chunks to make STAs from all possible contiguous subsets
% (of a specified minimum size).
% The results are stored in two variables, 'spatial_sensitivity' and 'combos'.
% combos{ss} is a vector of which chunks were used, and spatial_sensitivity{ss} is the 
% thresholded sum of STAs from those chunks.
%
% These two variables should be saved to a mat file 'pieces_analysis.mat'.
%
% The results can be viewed using the function temporal_subset_visualize.m
% The STAs can be combined together to form a normal datarun struct using temporal_subset_combine.m


% parameters
plot_figs = false;
start_combo = 1;


% choose source folders
clear sources
switch tsc
    case 1 % blueberry, 20 min chunks
        master = '2008-08-26-2/data001/data001';
        root_dir = '2008-08-26-2/data001/pieces/';
        for ss=1:8; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        %cell_spec = [47 378 1186 1923 2209 2358 ]; %{1,2,3,4,5,6};
        cell_spec = {1,2,3,4,5,6};
    case 2 % plantain, 15 min chunks
        master = '2008-08-27-5/data003/data003/data003';
        root_dir = '2008-08-27-5/data003/data003/pieces/';
        for ss=1:8; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        cell_spec = {1,2,3,4,5,6};
    case 3 % kiwi, 15 min chunks
        master = '2008-05-13-3/data006/data006';
        root_dir = '2008-05-13-3/data006/pieces/';
        for ss=1:9; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        cell_spec = 'all';
    case 4 % butterfly, 15 min chunks
        master = '2008-12-12-1/data005/data005';
        root_dir = '2008-12-12-1/data005/pieces/';
        for ss=1:8; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        cell_spec = 'all';
    case 5 % cherry, 10 min chunks
        master = '2009-02-28-0/data006/data006';
        root_dir = '2009-02-28-0/data006/pieces/';
        for ss=1:12; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 4;
        cell_spec = 'all';
    case 6 % pomegranate, 10 min chunks
        master = '2007-08-21-1/data003/data003';
        root_dir = '2007-08-21-1/data003/pieces/';
        for ss=1:6; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        cell_spec = 'all';
    case 7 % cherimoya, 10 min chunks
        master = '2008-03-25-3/data002/data002';
        root_dir = '2008-03-25-3/data002/pieces/';
        for ss=1:7; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 3;
        cell_spec = 'all';
        
    case 8 % mango, 10 min chunks
        master = '2008-04-30-2/data004/data004/data004';
        root_dir = '2008-04-30-2/data004/data004/pieces/';
        for ss=1:22; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 4;
        cell_spec = 'all';
        start_combo = 162;
    case 9 % grapes, 10 min chunks
        master = '2007-03-27-2/data014/data014/data014';
        root_dir = '2007-03-27-2/data014/data014/pieces/';
        for ss=1:23; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 4;
        cell_spec = 'all';
        start_combo = 165;
    case 10 % peach, 10 min chunks
        master = '2008-08-27-0/data001/data001';
        root_dir = '/lightning/Analysis/gauthier/2008-08-27-0/data001/pieces/';
        for ss=1:16; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 4;
        cell_spec = 'all';
        start_combo = 87;
    case 11 % apricot, 10 min chunks
        master = '2009-04-13-5/data005/data005';
        root_dir = '/lightning/Analysis/gauthier/2009-04-13-5/data005/pieces/';
        for ss=1:12; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 4;
        cell_spec = 'all';
    case 12 % plum, 15 min chunks
        master = '2008-04-22-5/data006/data006';
        root_dir = '2008-04-22-5/data006/pieces/';
        for ss=1:4; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = 2;
        cell_spec = 'all';
        
    case 0 % XXX, 10 min chunks
        master = '';
        root_dir = '/pieces/';
        for ss=1:XXXXXXX; sources{ss} = {sprintf('%spart%d/part%d',root_dir,ss,ss)}; end
        min_count = XXXXX;
        cell_spec = 'all';
end


% generate possible combinations
combos = cell(0);
for total_count = min_count:length(sources)
    for start_source = 1:(length(sources)-total_count+1)
        combos = {combos{:},[start_source:start_source+total_count-1]};
    end
end


% cycle through each possible combination
for cc = start_combo:length(combos)
    combo = combos{cc};
    tic

    % make a datarun struct with summed RFs
    
    % initialize with info from the master datarun
    datarun_sum = load_data(master);
    datarun_sum = load_sta(datarun_sum,'load_sta',[],'verbose',1);
    datarun_sum = load_params(datarun_sum,struct('verbose',1));
    
    % load basic info from each datarun
    datarun = load_data({sources{combo}});
    for dd=1:length(datarun)
        datarun{dd} = load_sta(datarun{dd},'load_sta',[],'verbose',1);
        datarun{dd} = load_neurons(datarun{dd});
    end

    % sum up STAs from each datarun, store in datarun_sum
    
    % identify number of frames
    num_frames = size(get_sta(datarun{1},datarun{1}.cell_ids(1)),4);
    
    % note which combo is being started
    fprintf('Starting combo %d of %d...\n',cc,length(combos))
    disp(combo)
    
    % go throuch each cell
    cell_indices = get_cell_indices(datarun_sum,cell_spec);
    for cell_index=cell_indices
        fprintf('.')
        cell_id = datarun_sum.cell_ids(cell_index);
        
        % initialize
        sta_temp = zeros(datarun_sum.stimulus.field_height,datarun_sum.stimulus.field_width,3,num_frames);
        % go through each datarun, add its version
        for dd =1:length(datarun)
            % check for cell id
            if any(datarun{dd}.cell_ids==cell_id)
                % get STA
                sta_partial = get_sta(datarun{dd},cell_id);
                % multiply by number of spikes
                sta_partial = sta_partial*length(datarun{dd}.spikes{get_cell_indices(datarun{dd},cell_id)});
                % add it in
                sta_temp = sta_temp + sta_partial;
            end
        end
        % convert to RF
        % get sig_stixels
        ss_params = struct;
        ss_params.strength = {'inner or',[0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};
        sig_stixels = significant_stixels(sta_temp,ss_params);
        % get RF
        rf = rf_from_sta(sta_temp,'sig_stixels',sig_stixels);
        % store in datarun_sum
        datarun_sum.stas.rfs{cell_index} = rf;
    end
    
    % compute spatial sensitivity
    
    % parameters
    clear spat_sens_params

    % how to combine RGB values of the RF
    spat_sens_params.strength = {'inner or',...
        [0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};

    % how to filter the RF before looking for significant stixels
    spat_sens_params.filter = [];

    % how to find significant stixels
    spat_sens_params.selection_params = struct('type','thresh','thresh',4.5);
    %params.selection_params = struct('type','max','thresh',5);

    % how to combine stixels between cells
    spat_sens_params.combine_stixels = 'sum';
    %cone_loc_params.combine_stixels = 'max';

    % online readout of what's going on
    spat_sens_params.verbose = true;
    spat_sens_params.fig_single_cell = [];
    spat_sens_params.foa_spat_sens = 102;
    
    % finally, compute it
    spatial_sensitivity{cc} = compute_spatial_sensitivity(datarun_sum, cell_spec, spat_sens_params);
    
    % plot it
    if plot_figs
        figure(cc);clf;imagesc(spatial_sensitivity{cc});axis image;drawnow
    end
        
    fprintf('\n   finished in %0.1f min...  ',round(toc)/60)
    
    % save somewhere
    save([server_path root_dir 'pieces_analysis' num2str(cc)],'spatial_sensitivity','combos')

    fprintf('saved.\n')
end

% view
for ss=1:length(spatial_sensitivity);figure(1);clf;imagesc(spatial_sensitivity{ss});disp(combos{ss});pause;end

