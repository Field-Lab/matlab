% generator signal --nishal


date='2015-04-14-2';
concatname='data000';
stim_time = 1800; % seconds
% Wrong Movie Information
file_name = [date, '/', concatname, '/', concatname];
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111-40x40.xml';

cell_type = {'ON parasol'};
num_frames = 30; % both have to be run with the name number of frames

cell_specification = 2538; %ON parasol
vision_id = 2538;
%% END OF INPU

% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames =1:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_types{1} = 'ON parasol';
% NOTE: cell_types{i}.name to be chosen to include the target cell of interest. 
addpath(genpath('../../../code/'));
datarun = get_sta_summaries(datarun, datarun.cell_types{1}.name, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
    
movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111-40x40.xml' %RGB-8-1-0.48-11111.xml';%RGB-8-1-0.48-11111-80x40
datarun = load_java_movie(datarun, movie_spec);
%datarun = get_snls(datarun, datarun.cell_ids(get_cell_indices(datarun, datarun.cell_types{1}.name)),'frames',-18:0,'stimuli',[],'new',true);
datarun = get_snls(datarun,vision_id,'frames',-18:0,'stimuli',[],'new',true);
% NOTE: TODO : have to Vary significant stixels, STA length (in frames parameter which is set at -18 right now !)and see how shape of non-linearity changes!

[cell_numbers] = get_cell_indices( datarun, vision_id );

cellID=cell_numbers;
gen=datarun.stas.snls{cellID}.gen_signal;
spks=datarun.stas.snls{cellID}.spikes;

figure;
subplot(2,2,1);
hist(spks);
subplot(2,2,2);
scatter(gen,spks);
hold on;
x=[-1:0.01:1];
N=@(x) exp(datarun.stas.snls{cellID}.fit_params.a*x +datarun.stas.snls{cellID}.fit_params.b);
plot(x,N(x),'r');


%% added
spikes = full(spks);
data = [gen, spikes];

data_sorted = sortrows(data);
range = data_sorted(end,1) - data_sorted(1,1);

num_bins = 40;

%% bins even in number of spikes
spikes_per_bin = sum(spikes)/num_bins;

spike_count = 0;
bin_counter = 1;
start_iter = 1;
while bin_counter < num_bins
    for i = start_iter:size(data_sorted,1)
        spike_count = spike_count + data_sorted(i,2);
        if spike_count >= spikes_per_bin

            bin_edges(bin_counter) = data_sorted(i,1);
            total_spike(bin_counter) = spike_count;
            bin_counter = bin_counter + 1;
            spike_count= 0;
            break
        end
    end
    start_iter = i+1
end


bin_edges(num_bins) = data_sorted(end,1);
total_spike(num_bins) = sum(spikes) - sum(total_spike);

bin_edges_with_first = [data_sorted(1,1), bin_edges];

bin_centers = bin_edges_with_first(1:end-1) + (bin_edges_with_first(2:end) - bin_edges_with_first(1:end-1))/2;


width_bins = bin_edges_with_first(2:end) - bin_edges_with_first(1:end-1);
FR = total_spike./width_bins;

figure;
plot(bin_centers,FR, 'o-')



%% evenly spaced bins
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