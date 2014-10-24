function datarun = cone_preprocess(piece, data, cones, stim)

%% set path info
% i need to know:
% the name of the folder with the cone mappings
% the name of the stim file

% example of inputs:
%piece = '2012-01-27-1';
%data = 'data004';
%cones = '_snle_acquisition_2012-01-27-1_data004_data004-bayes-msf_15.00-RGB-2-5';
%stim = 'BW-2-5-0.48-22222-300x300-60.35.xml';

dataPath = strcat('/marte/snle/acquisition/',piece,'/',data,'/',data);
conePath = strcat('/marte/snle/lab/Experiments/Array/Shared/one/',cones);
stimPath = strcat('/marte/snle/acquisition/movie-xml2/',stim);
savePath = strcat('/Users/Vision/Desktop/freeman/data/',piece,'/',data,'/');

%% load data from java files into matlab
datarun = load_data(dataPath);
datarun.names.movie_xml_path = stimPath;
datarun = load_neurons(datarun);
datarun = load_java_movie(datarun);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun,'load_sta', 'all','save_sta',false);

if all(cellfun(@isempty, datarun.stas.rfs))
    datarun = get_sta_summaries(datarun, {1,2,3,4,5}, 'marks_params', struct('robust_std_method', 6));
end


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
fprintf('Computing cone input in frames %d to %d...', start_stim, end_stim)
T = text_waitbar;

field_width = datarun.stimulus.field_width;
field_height = datarun.stimulus.field_height;

% Pre-permute Wc to avoid permuting every frame...
ind = 1:size(Wc,1);
ind = reshape(ind, field_height, field_width, 3);
ind = permute(ind, [3 2 1]);
ind = reshape(ind, [], 1);
Wcp = Wc(ind,:);

% cycle through each stimulus
cone_inputs = zeros(length(stims), size(Wc,2), 'single');
for ss = 1:length(stims)
    T = text_waitbar(T, ss/length(stims) - 0.01);

    % note which stim
    this_stim = stims(ss);

    % get new frame
    STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
    new_frame = double(STAFrame.getBuffer());
    
    % convert to cone space
    cone_inputs(ss,:) = single((Wcp' * new_frame) - 0.5);
end
datarun.cone_inputs = cone_inputs;


%% OLD reexpress movie in terms of cone signals
% fprintf('Computing cone input in frames %d to %d...',start_stim,end_stim)
% T = text_waitbar;
% 
% % initialize storage variables
% field_width = datarun.stimulus.field_width;
% field_height = datarun.stimulus.field_height;
% 
% % cycle through each stimulus
% cone_inputs = zeros(length(stims), size(Wc,2), 'single');
% for ss = 1:length(stims)
%     T = text_waitbar(T, ss/length(stims) - 0.01);
% 
%     % note which stim
%     this_stim = stims(ss);
% 
%     % get new frame
%     STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
%     new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;
%     new_frame = reshape(new_frame, [], 1);
% 
%     % convert to cone space
%     cone_inputs(ss,:) = single(full(Wc' * double(new_frame)));
% end
% datarun.cone_inputs = cone_inputs;


%% get the spikes in each bin
datarun = get_binned_spikes(datarun, start_stim, end_stim, stims);

%% save out the file
%save(fullfile(savePath,'datarun.mat'), 'datarun', '-v7.3');