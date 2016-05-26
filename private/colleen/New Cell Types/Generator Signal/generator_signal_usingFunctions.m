% produce generator signal using functions written by Peter Li
clear
clc
close all


%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% date='2005-04-26-1';
% concatname='data006';
% file_name = [date, '/', concatname, '/',  concatname];
% movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
%
% cell_type = {'OFF-parasol'};
% num_frames = 30;
% cell_specification = 182; % visionID
% num_bins = 20; % Number of bins to use for the generator vs spikes graph


date='2016-02-17-1';
concatname='data000';
file_name = [date, '/', concatname, '/',  concatname];
movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-119.5.xml';

cell_type = {'ON parasol'};
num_frames = 30;
cell_specification = 1456; % visionID
num_bins = 10; % Number of bins to use for the generator vs spikes graph

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath= ['/Users/colleen/Desktop/Generator Signal/', date, '/', concatname, '/'];


% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames =1:30;% have to input as a vector list of frames, not the number of frames total
datarun=load_data(datarun,opt);


% Changed the threshold from 5 to 4 because wasn't finding any stixels
datarun = get_sta_summaries(datarun, cell_specification, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
    'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
    'thresh',4,'robust_std_method',1));

datarun = load_java_movie(datarun, movie_spec);

% Get the SNLs based on the significant stixels and the last 19 frames
datarun = get_snls(datarun,cell_specification,'frames',-29:0,'stimuli',[],'new',true);

[cellID] = get_cell_indices( datarun, cell_specification);

gen=datarun.stas.snls{cellID}.gen_signal;
spks=datarun.stas.snls{cellID}.spikes;

% Look at the spike pattern and the generator signal
% figure;
% subplot(2,
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



%% bins even in number of spikes
% spikes_per_bin = sum(spikes)/num_bins;
%
% spike_count = 0;
% bin_counter = 1;
% start_iter = 1;
% included_gen_signals = zeros(1,num_bins);
% total_spike = zeros(1,num_bins);
% while bin_counter < num_bins
%     for i = start_iter:size(data_sorted,1)
%         spike_count = spike_count + data_sorted(i,2);
%         if spike_count >= spikes_per_bin
%             included_gen_signals(bin_counter) = i; %Used to calculate bin mean
%             total_spike(bin_counter) = spike_count; % How many spikes are in each bin (~= spikes_per_bin)
%             bin_counter = bin_counter + 1;
%             spike_count= 0;
%             break
%         end
%     end
%     start_iter = i+1;
% end
%
% % Fill in the last element of total_spike with the remaining spikes (could be slightly less than spikes_per_bin
% total_spike(num_bins) = sum(spikes) - sum(total_spike);
%
% % In order to plot the binned generator signal, determine the mean of the
% % generator signal for each bin as well as the range of generator signals
% % in a bin
% bin_mean(1) = mean(data_sorted(1:included_gen_signals(1)));
% width_bins(1) = data_sorted(included_gen_signals(1),1) - data_sorted(1,1);
% for i = 2:num_bins-1
%     bin_mean(i) = mean(data_sorted(included_gen_signals(i-1):included_gen_signals(i)));
%     width_bins(i) = data_sorted(included_gen_signals(i),1) - data_sorted(included_gen_signals(i-1),1);
%
% end
% bin_mean(num_bins) = mean(data_sorted(included_gen_signals(num_bins-1)+1:size(data_sorted,1)));
% width_bins(num_bins) = data_sorted(end,1) - data_sorted(included_gen_signals(num_bins-1),1);
%
% % Normalize the spike count per bin by the width of the bin
% FR = total_spike./width_bins;
%
% figure;
% plot(bin_mean,FR, 'o-')



% %% evenly spaced bins
% bin_edges = linspace(data_sorted(1,1), data_sorted(end,1), num_bins+1);
% [P] = histcounts(data_sorted(:,1), bin_edges)
% count = 0;
% FR = zeros(1, num_bins);
% for i = 1:length(P)-1
%     start = P(i);
%     FR(i) = sum(data_sorted((count+1):(count+start), 2));
%     count = count + start;
% end
% bin_centers = bin_edges(1:end-1) + (bin_edges(2:end) - bin_edges(1:end-1))/2;
%
% figure; plot(bin_centers, FR)


%% Save number of generator signal values

num_of_gen_signals = size(data_sorted,1);
gen_signals_per_bin = num_of_gen_signals/num_bins;


bin_edges = [1:ceil(gen_signals_per_bin):size(data_sorted,1), size(data_sorted,1)];
binned_gen_signals=nan(num_bins,1);
mean_spikes=nan(num_bins,1);

for i = 1:length(bin_edges)-1
    binned_gen_signals(i) =data_sorted(bin_edges(i)) % mean(data_sorted(bin_edges(i:i+1),1));
    %     binned_spikes(i) = sum(data_sorted(bin_edges(i):bin_edges(i+1),2))
    mean_spikes(i) = mean(data_sorted(bin_edges(i):bin_edges(i+1),2));
end
figure;
plot(binned_gen_signals, mean_spikes, 'o-')
xlabel('Generator Signal')
ylabel('Spike Rate (spikes/bin)')
title({[date, ' ', concatname]; cell_type{1}; ['Cell ' num2str(cell_specification)]})


% Save the Resultf ~exist(filepath)
mkdir(filepath)

export_fig([filepath, 'Cell_',num2str(cell_specification)], '-pdf')
