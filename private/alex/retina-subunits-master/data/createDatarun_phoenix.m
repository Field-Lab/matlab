function createDatarun_pheonix(date,data,cones,stim)

%% set path info
% i need to know:
% the name of the folder with the cone mappings
% the name of the stim file

% example of inputs:
%date = '2012-01-27-1';
%data = 'data004';
%cones = '_snle_acquisition_2012-01-27-1_data004_data004-bayes-msf_15.00-RGB-2-5';
%stim = 'BW-2-5-0.48-22222-300x300-60.35.xml';

dataPath = strcat('/marte/snle/acquisition/',date,'/',data,'/',data);
conePath = strcat('/marte/snle/lab/Experiments/Array/Shared/one/',cones);
stimPath = strcat('/marte/snle/acquisition/movie-xml2/',stim);
savePath = strcat('/Users/Vision/Desktop/freeman/data/',date,'/',data,'/');

%% load data from java files into matlab
datarun = load_data(dataPath);
datarun.names.movie_xml_path = stimPath;
datarun = load_neurons(datarun);
datarun = load_java_movie(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun,'load_sta', 'all','save_sta',false);
datarun = get_sta_summaries(datarun, {1,2,3,4,5});

%% get ready to load the movie
start_time = 0;
end_time = datarun.duration; % change to load a subset...
time_offset = 0;
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
stims = start_stim:end_stim;

%% load cone info
datarun = import_single_cone_data(datarun,conePath);
load(fullfile(conePath,'Wc.mat'));

%% reexpress movie in terms of cone signals
fprintf('Computing cone input in frames %d to %d...',start_stim,end_stim)
T=text_waitbar;
start_time_ = clock; % note when it started
% initialize storage variables
field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;
cone_inputs = zeros(length(stims),size(Wc,2));
% cycle through each stimulus
for ss = 1:length(stims)
    T=text_waitbar(T,ss/length(stims) - 0.01);
    % note which stim
    this_stim = stims(ss);
    % get new frame
    STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
    new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;
    new_frame = reshape(new_frame,[],1);
    % convert to cone space
    cone_inputs(ss,:) = full(Wc'*double(new_frame));
end
datarun.cone_inputs = cone_inputs;

%% get the spikes in each bin
datarun = getBinnedSpikes(datarun,start_stim,end_stim,stims);

%% save out the file
save(fullfile(savePath,'datarun.mat'),'datarun','-v7.3');