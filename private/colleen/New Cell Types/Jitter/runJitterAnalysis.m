% runJitterAnalysis

clear
close all

dataparam.date='2015-09-23-7';
dataparam.concatname='d19-39/data029-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-2-0.48-11111.xml';
fitparam.num_frames = 30;

% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/data029-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039'];

% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
    dataparam.cell_specification = [530] %ON parasol
end
dataparam.cell_type = {'ON parasol'};
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

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);




%% Cell indicies

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
    dataparam.cell_specification = datarun.cell_types{cell_type_index}.cell_ids;
end


triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(dataparam.mdf_file,...
    triggers, 1,2);
%
[mvi] = load_movie(dataparam.mdf_file, triggers);
cell_indices = get_cell_indices(datarun, dataparam.cell_specification);
cell_ids=datarun.cell_ids(cell_indices);
for j = 1:length(cell_indices)
    stixel_size = 20;
        stixel_width = 20;
            stixel_height = 20;
    num_frames = 15;
    spikes=  datarun.spikes{cell_indices(j)};
    seed = 11111;
    state = Init_RNG_JavaStyle(seed);
    jitterX = nan(duration,1);
    jitterY = nan(duration,1);
    
    for i = 1:duration
        jitterX(i) = [mod(double(random_uint16(state)), stixel_width) - stixel_width/2];
        jitterY(i) = [mod(double(random_uint16(state)), stixel_height) - stixel_height/2];
    end
    save ('jitterX', 'jitterX');
     save ('jitterY', 'jitterY');
%     jitterX = zeros(size(jitterX));
%     jitterY = zeros(size(jitterY));
    num_colors =3;
    [sta] = compute_jitter_sta(datarun, dataparam.mdf_file, num_frames, spikes, jitterX, jitterY, stixel_size, num_colors);
    
    save(['/Users/colleen/Desktop/Jitter/2015-09-23-7/data028/Cell ', num2str(cell_ids(j))], 'sta')
end
