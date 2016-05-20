% ruvnJitterAnalysis
% git hash 75ce3d76b7b20220cd94d83eb96e8e3ec90dd20e

clear
close all
profile on
% dbstop if error
dataparam.date='2016-04-21-1';
dataparam.concatname='data005';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-40-1-0.48-11111-20x15-60.35.xml';
dataparam.stixel_size = 40;
dataparam.interval = 1;
dataparam.seed = 11111;
dataparam.refresh_rate = 60.35;
fitparam.num_frames = 30;
dataparam.x_dim = 800;
dataparam.y_dim = 600;
num_colors =3;
dataparam.num_of_interval = 50; % number of * in progress bar
frame_width = dataparam.x_dim/dataparam.stixel_size;
frame_height = dataparam.y_dim/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;



% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];

% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
  dataparam.cell_specification = [35] %ON parasol

end
dataparam.cell_type = {'all'};
%% END OF INPUT
% dataparam.folder = dataparam.cell_type{1};
% file path to save data and pictures
% dataparam.filepath=['/Users/colleen/Desktop/Fitting/',dataparam.date,'/',dataparam.concatname,'/data023/'];
% if ~exist([dataparam.filepath,dataparam.folder],'dir')
%     mkdir([dataparam.filepath,dataparam.folder]);
% end

% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];



%% Load Data2
slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);




%% Cell indicies
tic
% Find the type of the inputted cells
cell_type_index= zeros(1,size(dataparam.cell_type,2));
for num_cell_types = 1:size(dataparam.cell_type,2)
    for i = 1:size(datarun.cell_types,2)
        right_cell_type = strcmpi(datarun.cell_types{i}.name, dataparam.cell_type{num_cell_types}); % case insensitive
        if right_cell_type == 1;
            cell_type_index(num_cell_types) = i;
            break
        end
        cell_type_index(num_cell_types) = 0;% couldn't find the right cell type
    end
    
end

% Set the cell_specification to all the cell of the inputted type
if select_cells ~= 1
    dataparam.cell_specification = datarun.cell_ids;%datarun.cell_types{cell_type_index}.cell_ids;
end


% triggers=datarun.triggers(1:500); %onsets of the stimulus presentation

triggers = [datarun.triggers; [datarun.triggers(end) + mean(diff(datarun.triggers)):mean(diff(datarun.triggers)):datarun.triggers(end) + 600*mean(diff(datarun.triggers))]'];
[~,height,width,duration,refresh] = get_movie_ath(dataparam.mdf_file,triggers, 1,2);
%
[mvi] = load_movie(dataparam.mdf_file, triggers);
cell_indices = get_cell_indices(datarun, dataparam.cell_specification);
cell_ids=datarun.cell_ids(cell_indices);
stixel_size = dataparam.stixel_size;
stixel_width = dataparam.stixel_size;
stixel_height = dataparam.stixel_size;
seed = dataparam.seed;

for j = 1:length(cell_indices)
    
    
    spikes{j}=  datarun.spikes{cell_indices(j)};
    
end

% jitterX = load('jitterX_16.mat');
% jitterX = jitterX.jitterX;
% jitterY = load('jitterY_16.mat');
% jitterY = jitterY.jitterY;

if exist(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname, '/jitterX.mat']) ~=0 && exist(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname, '/jitterY.mat']) ~=0
    
    jitterX = load(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname,'/', 'jitterX.mat']);
    jitterX = int16(jitterX.jitterX);
    jitterY = load(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname,'/', 'jitterY.mat']);
    jitterY = int16(jitterY.jitterY);
else
    
    
    
    
    state = Init_RNG_JavaStyle(seed);
    jitterX = int16(zeros(duration,1));
    jitterY = int16(zeros(duration,1));
    
    
    for i = 1:duration
        jitterX(i) = [mod(double(random_uint16(state)), stixel_width) - stixel_width/2];
        jitterY(i) = [mod(double(random_uint16(state)), stixel_height) - stixel_height/2];
    end
    
    
    if ~exist(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname, '/'])
        mkdir(['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname,'/'])
    end
    
    save (['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname,'/', 'jitterX'], 'jitterX');
    save (['/Volumes/Lab/Users/crhoades/Jitter_Shifts/', dataparam.date, '/', dataparam.concatname,'/', 'jitterY'], 'jitterY');
end





[sta] = compute_jitter_sta_opt(datarun, dataparam.mdf_file, fitparam.num_frames, spikes, jitterX, jitterY, stixel_size, num_colors, dataparam);
for i = 1:size(spikes,2)
    temp = squeeze(sta(:,:,:,:,i));
    if ~exist(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/', dataparam.concatname])
        mkdir(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/',dataparam.concatname]);
    end
    save(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/', dataparam.concatname, '/Cell ', num2str(cell_ids(i))], 'temp')
end

toc
profile off
profile viewer
[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);
