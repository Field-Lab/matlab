function [dataset,stds_before] = spike_plot_and_select_spikes(dataset,projection_type,electrodes_for_projections,figures_to_use,proj_index)

% perform PCA on a subset of spikes
dataset.electrodes_for_pca = electrodes_for_projections;
dimension_prefix = 'PC';

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

%     DISCARDING ELECCTRODES
    [spikes_weighted,junk,spike_waveform_points_to_use] = ...
        spike_weight_calculator(spikes_original,'electrodes',electrodes_for_projections,length(dataset.electrodes),dataset.window_length);



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
    %proj_struct.selected_indices = spike_select_spikes(proj_struct,1,'b');
    
    % note which spikes were selected
    %selected_spike_waveforms = proj_struct.spike_waveforms(proj_struct.selected_indices{1},:);
    %selected_spike_times = proj_struct.spike_times(proj_struct.selected_indices{1},1);

    % plot statistics about these spikes
    %spike_compute_and_plot_acf(selected_spike_times,acf_range,proj_struct.plot_axes.acf_axis{1});
    %spike_plot_average_spike(selected_spike_waveforms,length(electrodes_for_projections),dataset.window_length,proj_struct.plot_axes.avg_spike_axis{1});
    %spike_compute_and_plot_sta(selected_spike_times,dataset.triggers,dataset.mdf_file,proj_struct.plot_axes.sta_axis{1});
    
    
    % save in dataset variable
    dataset.(projection_type){projection_index_number} = proj_struct;
    
end


