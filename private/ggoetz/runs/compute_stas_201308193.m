clear;
addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/private/ggoetz'));
datarunpath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data001/data001';
moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/eye-movement/eye-long-v2' ;

N_SPIKES_STA = 5000;

%% Load the data run

datarun = load_data(datarunpath);
datarun = load_neurons(datarun);

% RRS works in samples, but it is more natural to work in samples for us.
% Let's convert the datarun to samples.
datarun = convert_datarun_times_to_samples(datarun);

%% Get the time of the image refreshes from the ttls

t_frames = time_imrefresh_from_ttls(datarun.triggers);

%% Calculate the order in which the frames of the natural movie were shown

repeats = 60;
frames_a = 3600;
frames_b = 3600*2;
fo = frames_order_natural_movie(frames_a, frames_b, repeats);

%% Get a single cell STA

st = datarun.spikes{1};
st = st(1:N_SPIKES_STA);
staind = compute_ta_ind(st, t_frames);
staindremapped = remap_event_indices(staind, 1:length(fo), fo);
sta = compute_ta_from_ind(staindremapped, moviechunksfolder);