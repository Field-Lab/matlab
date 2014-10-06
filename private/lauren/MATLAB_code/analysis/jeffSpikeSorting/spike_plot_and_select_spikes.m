function [dataset,stds_before] = spike_plot_and_select_spikes(dataset,projection_type,electrodes_for_projections,figures_to_use,proj_index)

% perform PCA on a subset of spikes

% save parameters
switch projection_type
    case 'pca'
        dataset.electrodes_for_pca = electrodes_for_projections;
        dimension_prefix = 'PC';
    case {'lda', 'nnwpca'}
        dataset.electrodes_for_lda = electrodes_for_projections;
        dimension_prefix = 'LD';
    case 'wavelet'
        dataset.electrodes_for_pca = electrodes_for_projections;
        dimension_prefix = 'wavelet decomp.';
end

% SET DISPLAY PARAMETERS
bins = 230 ; %for density gradient
%dimensions_to_plot=[1 2;1 3;2 3]; %list of which plots to make
dimensions_to_plot=[1 2;1 3;2 3]; %list of which plots to make
acf_range = 30; %in msec

if isempty(electrodes_for_projections)
    electrodes_for_projections = 1:length(dataset.electrodes);
end

% PUT SPIKES IN NEW VARIABLE
spikes_original = dataset.raw_spikes;

% dynamic variance reduction
switch 0
    case 1
        constant = 1000;
        spikes_original = sqrt(abs(spikes_original)+constant).*sign(spikes_original);
    case 2
        constant = 50;
        change = abs(spikes_original)>constant;
        spikes_original(change) = (sqrt(abs(spikes_original(change))-constant)+constant).*sign(spikes_original(change));
    case 3
        constant = 5;
        spikes_original = sqrt(spikes_original.^2+constant).*sign(spikes_original);
    case 4  % variance reduction each window separately
        
        %make sure to skip the center electrode -- we do not want to reduce
        %variance there!
        reshaped_spikes = reshape(spikes_original(:,(2+dataset.window_length):end)',dataset.window_length,[]);
        stds_before = std(reshaped_spikes,[],1);
        switch 1
            case 1
                std_thresh = 40;
                norm_factors = stds_before;
                norm_factors(stds_before < std_thresh) = 1;
            case 2
        end
        reshaped_normed_spikes = reshaped_spikes./repmat(norm_factors,dataset.window_length,1);
        stds_after = std(reshaped_normed_spikes,[],1);
        figure(300);hist(stds_before,30);
        figure(200);plot(stds_before,stds_after+10*randn(1,length(stds_after)),'.')
        spikes_original(:,(2+dataset.window_length):end) = reshape(reshaped_normed_spikes,6*dataset.window_length,[])';
    case 4.5  % variance reduction each window separately (exclude center)
        reshaped_spikes = reshape(spikes_original(:,(2):end)',dataset.window_length,[]);
        stds_before = std(reshaped_spikes,[],1);
        switch 1
            case 1
                std_thresh = 40;
                norm_factors = stds_before;
                norm_factors(stds_before < std_thresh) = 1;
            case 2
        end
        reshaped_normed_spikes = reshaped_spikes./repmat(norm_factors,dataset.window_length,1);
        stds_after = std(reshaped_normed_spikes,[],1);
        figure(300);hist(stds_before,40);
        figure(200);plot(stds_before,stds_after+10*randn(1,length(stds_after)),'.')
        spikes_original(:,(2):end) = reshape(reshaped_normed_spikes,7*dataset.window_length,[])';
    case 5  % variance reduction of whole spike (all windows) but exclude center
        reshaped_spikes = reshape(spikes_original(:,(2+dataset.window_length):end)',dataset.samples_per_spike-dataset.window_length,[]);
        stds_before = std(reshaped_spikes,[],1);
        switch 1
            case 1
                std_thresh = 75;
                norm_factors = stds_before;
                norm_factors(stds_before < std_thresh) = 1;
            case 2
        end
        reshaped_normed_spikes = reshaped_spikes./repmat(norm_factors,dataset.samples_per_spike-dataset.window_length,1);
        stds_after = std(reshaped_normed_spikes,[],1);
        figure(200);plot(stds_before,stds_after+10*randn(1,length(stds_after)),'.')
        spikes_original(:,(2+dataset.window_length):end) = reshape(reshaped_normed_spikes,6*dataset.window_length,[])';
    case 6 % variance reduction of whole spike (all windows) but do not exclude center, threshold std
        reshaped_spikes = reshape(spikes_original(:,(2):end)',dataset.samples_per_spike,[]);
        stds_before = std(reshaped_spikes,[],1);
        switch 1
            case 1
                std_thresh = 50;
                norm_factors = stds_before;
                norm_factors(stds_before < std_thresh) = 1;
        end
        reshaped_normed_spikes = reshaped_spikes./repmat(norm_factors,dataset.samples_per_spike,1);
        stds_after = std(reshaped_normed_spikes,[],1);
        
        figure(200);plot(stds_before,stds_after+10*randn(1,length(stds_after)),'.')
        spikes_original(:,(2):end) = reshape(reshaped_normed_spikes,7*dataset.window_length,[])';
    case 7 % variance reduction of whole spike (all windows) but do not exclude center, based on cumulative variance histogram
        reshaped_spikes = reshape(spikes_original(:,(2):end)',dataset.samples_per_spike,[]);
        stds = std(reshaped_spikes,[],1);
        
        figure(666);
        switch 2
            case 1
                range = min(stds):0.1:max(stds);
                n_elements = histc(stds,range);
                c_elements = cumsum(n_elements);
                bar(range,c_elements);
            case 2
                plot(sort(stds),1:length(stds))
        end
        

        %compute scaling coefficients
        sorted = sort(stds);
        scaling_factors = zeros(length(stds),1);
         
        for ii = 1:length(stds)
           currIndex = find(sorted == stds(ii));
           scaling_factors(ii) = sum(stds(1:currIndex));
        end
         
        %obtain the variance
        scaling_factors = sqrt(sqrt(scaling_factors'));

        %before normalizing by the values obtained from the cumulative
        %variance histogram, we must normalize the variance of each spike
        %to one.
        norm_factors_primary = stds;
        reshaped_spikes = reshaped_spikes./repmat(norm_factors_primary,dataset.samples_per_spike,1);
         
         
        %now multiply by scaling coefficients
        reshaped_normed_spikes = reshaped_spikes.*repmat(scaling_factors,dataset.samples_per_spike,1);
        stds_after = std(reshaped_normed_spikes,[],1);
        
        figure(200);plot(stds,stds_after+0*randn(1,length(stds_after)),'.')
        spikes_original(:,(2):end) = reshape(reshaped_normed_spikes,7*dataset.window_length,[])';
end


% switch 1
%     case 1 % use entropy as the weighting function
%         spike_weight_mi = dataset.weights;
% end


% SELECT WHICH PART OF THE SPIKE WAVEEFORM TO USE BY EITHER:
switch 1
    
    case 1 %     DISCARDING ELECCTRODES
        [spikes_weighted,junk,spike_waveform_points_to_use] = ...
            spike_weight_calculator(spikes_original,'electrodes',electrodes_for_projections,length(dataset.electrodes),dataset.window_length);
        %   OR:
        
    case 2
        %     USING A WEIGHTING FUNCTION
        
        if 1  % vector of 1s and 0s
            spike_weight_mi = zeros(1,dataset.samples_per_spike);
            % choose which timepoints to use
            switch 2
                case 1 % set variance threshold
                    tp = var(dataset.raw_spikes(1:1000,2:end)) > 10^2.7;
                case 2 % just cetain timepoints
                    time_points = (9:13)';elecs = [1 2 4 6];
                    samples_per_elec = dataset.window_length;
                    tp = repmat(time_points,1,length(elecs))+repmat((elecs-1)*samples_per_elec + 1,length(time_points),1) - 1;
                case 3 % specially picked
                    tp = [7     9    11    12    13    14    15    16    40    41    42    43    44    45    46    49    50  51    52   110   112   113   114   144   145   146   147   148   149   151   152   153   154   155   156   177   178   180   181   182   185   213   214   215   219   220];
                case 4 % use the weighting function
                    tp = dataset.weights > .1;
            end
            % set those values to 1
            spike_weight_mi(tp) = 1;
            
        else  % vector of scalars
            spike_weight_mi = dataset.weights;
        end
            
        [spikes_weighted,junk,spike_waveform_points_to_use] = ...
            spike_weight_calculator(spikes_original,'samples',spike_weight_mi,length(dataset.electrodes),dataset.window_length);
end


% SHOW AVERAGE WAVEFORM OF ALL SPIKES, HIGHLIGHTING THOSE USED WITH RED DOTS
figure(2); title('average waveform of all spikes, those used highlighted with red ');
plot(spike_waveform_points_to_use-1,mean(spikes_original(:,spike_waveform_points_to_use),1),'r.');
hold on;
spike_plot_average_spike(spikes_original(:,2:size(spikes_original,2)), ...
    length(dataset.electrodes),dataset.window_length,gca);
hold off;


% COMPUTE AND PLOT PC PROJECTION AND SELECT SPIKES

% information is saved in fields of a struct array
% for example, for pca, these are the fields:
%
% dataset.pca{1}.spike_times
%               .spike_waveforms
%               .spike_projections
%               .selected_indices
%               .plot_axes.projection_axes{k} for k = 1 to 3
%                         .avg_spike_axis
%                         .acf_axis
%                         .sta_axis
%               .density_histogram_bins
%
% each subset gets a new number, so there will be
% dataset.pca{1}, dataset.pca{2}, dataset.pca{3}, ...
% it's up to the user to keep track of which is which and clear out old ones
%
% for lda, there would be another set like this:
% dataset.lda{1}, dataset.lda{2}, ...
%
% for information on what each of the fields is for, see below


for subset_number = 1:length(figures_to_use)
    
    % specify which element in the struct array this will be.
    % it will be used like this: dataset.pca(projection_index_number).spike_times
    projection_index_number = proj_index + subset_number - 1;
    
    % load parameters from dataset
    proj_struct.electrodes_for_projections = electrodes_for_projections;
    proj_struct.window_length = dataset.window_length;
    proj_struct.triggers = dataset.triggers;
    proj_struct.mdf_file = dataset.mdf_file;
    
    % all information will be stored in a temporary struct called proj_struct.  once
    % the information is entered, it will be saved in the dataset variable (see last step below)
    
    % set up axes in the figure
    proj_struct.plot_axes = spike_set_up_axes(figures_to_use(subset_number),1);
  
    % save the spikes_waveforms in the projection information
    switch subset_number
        case 1
            % the first time, use all spikes
            %proj_struct.spike_waveforms = spikes_original(:,2:size(spikes_original,2));
            proj_struct.spike_waveforms = spikes_original(:,2:end);
            proj_struct.spike_waveform_points_to_use = spike_waveform_points_to_use;
            proj_struct.spike_times = spikes_original(:,1);
        otherwise
            % in subsequent subsets, use only the spikes which were selected previously
            % first get all spikes from last time and which indices were selected
            spikes_from_last_time = dataset.(projection_type){projection_index_number - 1}.spike_waveforms;
            spike_times_from_last_time = dataset.(projection_type){projection_index_number - 1}.spike_times;
            indices_selected_last_time = dataset.(projection_type){projection_index_number - 1}.selected_indices{1};
            % then extract just the spikes which were selected
            proj_struct.spike_waveforms = spikes_from_last_time(indices_selected_last_time,:);
            proj_struct.spike_times = spike_times_from_last_time(indices_selected_last_time,1);
    end
    
    %compute projections and save in proj_struct
    proj_struct = spike_compute_projection(projection_type,proj_struct,dataset);
    
    %load parameters to plot the projections
    proj_struct.dimension_prefix = dimension_prefix;
    proj_struct.dimensions_to_plot = dimensions_to_plot;
    proj_struct.density_histogram_bins = bins;
    
    
    %plot projections in projection_axes
    proj_struct = spike_plot_density_gradient(proj_struct);

    % have user select spikes
    proj_struct.selected_indices = spike_select_spikes(proj_struct,1,'b');
    
    % note which spikes were selected
    selected_spike_waveforms = proj_struct.spike_waveforms(proj_struct.selected_indices{1},:);
    selected_spike_times = proj_struct.spike_times(proj_struct.selected_indices{1},1);

    % plot statistics about these spikes
    spike_compute_and_plot_acf(selected_spike_times,acf_range,proj_struct.plot_axes.acf_axis{1});
    spike_plot_average_spike(selected_spike_waveforms,length(electrodes_for_projections),dataset.window_length,proj_struct.plot_axes.avg_spike_axis{1});
    %spike_compute_and_plot_sta(selected_spike_times,dataset.triggers,dataset.mdf_file,proj_struct.plot_axes.sta_axis{1});
    
    
    % save in dataset variable
    dataset.(projection_type){projection_index_number} = proj_struct;
    
end


