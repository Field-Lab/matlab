% produce generator signal using functions written by Peter Li

function [plot_axes]  = computeAndPlot_genSignal(datarun, cell_specification, run_opts )

% Changed the threshold from 5 to 4 because wasn't finding any stixels
datarun = get_sta_summaries(datarun, cell_specification, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
    'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
    'thresh',4,'robust_std_method',1));

datarun = load_java_movie(datarun, run_opts.movie_spec);

% Get the SNLs based on the significant stixels and the last 19 frames
datarun = get_snls(datarun,cell_specification,'frames',-(run_opts.num_frames-1):0,'stimuli',[],'new',true);

[cellID] = get_cell_indices( datarun, cell_specification);
vision_cell_ids=datarun.cell_ids(cellID);
range =  -0.75:0.1:0.75;
y = nan(length(range), length(cellID));
for j = 1:length(cellID)
    gen=datarun.stas.snls{cellID(j)}.gen_signal;
    spks=datarun.stas.snls{cellID(j)}.spikes;
    
    % Look at the spike pattern and the generator signal
    % figure;
    % subplot(2,2,1);
    % hist(spks);
    % subplot(2,2,2);
    % scatter(gen,spks);
    % hold on;
    % x=[-1:0.01:1];
    % % Plot the exponential fit
    % N=@(x) exp(datarun.stas.snls{cellID}.fit_params.a*x +datarun.stas.snls{cellID}.fit_params.b);
    % plot(x,N(x),'r');
    
    
    %% Plot the generator signal against the spike count, binning similar generator signals together
    spikes = full(spks);
    data = [gen, spikes];
    % Sort the data from lowest generator signal to highest. Column 2 is the
    % corresponding spike count for each generator signal
    data_sorted = sortrows(data);
    
    %% Save number of generator signal values
    
    num_of_gen_signals = size(data_sorted,1);
    gen_signals_per_bin = num_of_gen_signals/run_opts.num_bins;
    
    
    bin_edges = [1:ceil(gen_signals_per_bin):size(data_sorted,1), size(data_sorted,1)];
    binned_gen_signals=nan(run_opts.num_bins,1);
    mean_spikes=nan(run_opts.num_bins,1);
    
    for i = 1:length(bin_edges)-1
        binned_gen_signals(i) =mean(data_sorted(bin_edges(i:i+1),1)); % Plot x axis from left edge of bin: data_sorted(bin_edges(i));
        mean_spikes(i) = mean(data_sorted(bin_edges(i):bin_edges(i+1),2));
    end
    
    normalized_mean_spikes = (mean_spikes - min(mean_spikes))/(max(mean_spikes)-min(mean_spikes));
    % figure;
    % h = plot(binned_gen_signals, mean_spikes, 'o-');
    % xlabel('Generator Signal')
    % ylabel('Spike Rate (spikes/bin)')
    % title({[run_opts.date, ' ', run_opts.concatname]; run_opts.cell_type{1}; ['Cell ' num2str(vision_cell_ids(j))]})
    
    BGS_interp= interp1(binned_gen_signals, normalized_mean_spikes,range);
    y(:, j) = BGS_interp;
    
end

params.foa = 0;
params.clear= 0;
plot_axes = set_up_fig_or_axes(params.foa,params.clear);

plot(range, mean(y,2), 'o-', 'Parent', plot_axes);
xlabel('Generator Signal')
ylabel('Normalized Average Spike Rate')
title({[run_opts.date, ' ', run_opts.concatname]; run_opts.cell_type{1}})


% Save the Result
% if ~exist(run_opts.filepath)
%     mkdir(run_opts.filepath)
% end
% export_fig([run_opts.filepath, 'Cell_',num2str(cell_specification)], '-pdf');
