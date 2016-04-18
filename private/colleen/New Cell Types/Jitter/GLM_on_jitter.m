%% GLM on Jitter Data

clear
close all

dbstop if error
dataparam.date='2016-02-17-6';
dataparam.concatname='data026';
dataparam.concatname_ref='data025';

dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-2-0.48-22222-119.5.xml';
fitparam.num_frames = 20;
frame_width = 640/16;
frame_height = 320/16;
stixels_per_frame = frame_width*frame_height;

num_colors =3;
% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_ref = [dataparam.date, '/', dataparam.concatname_ref,'/', dataparam.concatname_ref];

dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];
% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
    %        dataparam.cell_specification = [184] %ON parasol
    %     dataparam.cell_specification = [481 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 6484 185 1217 1221 1938 1942 2012 3037 3904 3908 4112 4893 4982 5523 5526 6094 6533 7056 7475 5061 3618 3811] %ON parasol
    dataparam.cell_specification = 1936;
    dataparam.cell_specification_ref = 1936;
    
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

datarun_ref.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_ref, '.neurons'];
datarun_ref.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_ref, '.params'];
datarun_ref.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_ref, '.sta'];





%% Load Data2
slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

datarun_ref=load_data(datarun_ref,opt);



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
cell_indices_ref = get_cell_indices(datarun_ref, dataparam.cell_specification_ref);

cell_ids=datarun.cell_ids(cell_indices);
cell_ids_ref=datarun.cell_ids(cell_indices_ref);

stixel_size = 16;
stixel_width = 16;
stixel_height = 16;
seed = 22222;

jitterX = load('jitterX_16.mat');
jitterX = jitterX.jitterX;
jitterY = load('jitterY_16.mat');
jitterY = jitterY.jitterY;
for j = 1:length(cell_indices)
    
    
    spikes{j}=  datarun.spikes{cell_indices(j)};
    %
    %     state = Init_RNG_JavaStyle(seed);
    %     jitterX = nan(duration,1);
    %     jitterY = nan(duration,1);
    %     duration
    %     for i = 1:duration
    %         if mod(i,1000) == 1
    %             disp(num2str(i));
    %         end
    %         for trash = 1:stixels_per_frame
    %             a = random_uint16(state);
    %         end
    %
    %         %         if mod(i, stixels_per_frame) == 0
    %         jitterX(i) = [mod(double(random_uint16(state)), stixel_width) - stixel_width/2];
    %         jitterY(i) = [mod(double(random_uint16(state)), stixel_height) - stixel_height/2];
    %
    %
    %
    %
    %     end
    %
    %     save ('jitterX', 'jitterX');
    %     save ('jitterY', 'jitterY');
    %
    %
    
    
    
end


% adjust for glitch
for i = 1:size(spikes,2)
    ind = find(spikes{i} > 92851/60 & spikes{i} <= 92900/60);
    if ~isempty(ind)
        spikes{i}(ind(end)+1:end) = spikes{i}(ind(end)+1:end) - (92900/60-92851/60);
        spikes{i} = [spikes{i}(1:ind(1)-1); spikes{i}(ind(end)+1:end)];
        
    end
    
    
end


% blocked spikes should be a cell of the spike times WITHIN each block
%   so each block's time should start over from 0

% Blocks are 200 frames at 120Hz = 1.66 seconds
cutoff = 1.66;
counter1 = 1;
counter2 = 1;
for i = 1:length(spikes{1,1})
    if spikes{1,1}(i) < 1.66
        blocked_spikes{counter1}(counter2,1) = spikes{1,1}(i);
        counter2 = counter2+1;
    else
        spikes{1,1} = spikes{1,1} - 1.66;
        counter1 = counter1+1;
        counter2 = 1;
    end
    
end

%% Get movie blocks




save_path = '/Volumes/Lab/Users/crhoades/JitterMovie/2016-02-17-6/data026/';
datarun.triggers = [datarun.triggers; [datarun.triggers(end) + mean(diff(datarun.triggers)):mean(diff(datarun.triggers)):datarun.triggers(end) + 300*mean(diff(datarun.triggers))]'];
%
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-2-0.48-22222-119.5.xml';

triggers = datarun.triggers;
% [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
%     triggers, 1,2);
bw = 0;
[inputs, refresh, duration] = get_wn_movie_ath(datarun, mdf_file, bw);
% test = reshape(1:800*85000, 800, 85000);

% test = [repmat([1;2;3],800,1), repmat([4;5;6],800,1),repmat([7;8;9],800,1),repmat([10;11;12],800,1)]
% inputs = test;
% real_frame(:,1) = inputs(1:(25));

% account for dropped frames
% data026
% 1857 is x coordinate of sharp peak in diff(triggers) graph
triggers = [triggers(1:1857); triggers(1858:end) - mean(diff(triggers(1:1857)))];
jitter_x_new = [jitterX(1:92851); jitterX(92900:end)]; %1857*100/2 (refresh is 2)+1 : 1857*100/2 (refresh is 2)+50
jitter_y_new = [jitterY(1:92851); jitterY(92900:end)];
inputs = [inputs(:,1:92851), inputs(:,92900:end)];

%% data024
% triggers = [triggers(1:392); triggers(393:end) - mean(diff(triggers(1:392)))];
% jitter_x_new = [jitter_x(1:19601); jitter_x(19650:end)]; %392*100/2 (refresh is 2)+1 : 392*100/2 (refresh is 2)+50
% jitter_y_new = [jitter_y(1:19601); jitter_y(19650:end)];
% inputs = [inputs(:,1:19601), inputs(:,19650:end)];
% for i = 1:size(spikes,2)
%     ind = find(spikes{i} > 19601/60 & spikes{i} <= 19650/60);
%     spikes{i}(ind(2)+1:end) = spikes{i}(ind(2)+1:end) - (19650/60-19601/60);
%     spikes{i} = [spikes{i}(1:ind(1)-1); spikes{i}(ind(2)+1:end)];
%
% end




image_width = 40;
image_height = 20;

real_frame = zeros(image_width, image_height, 3);
real_frame(:,:,1,1) = reshape(inputs(1:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,2,1) = reshape(inputs(2:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,3,1) = reshape(inputs(3:3:image_width*image_height*3)',image_width, image_height);

pointer = image_width*image_height*3+1;
%     pointer = 2+25+2;
i =2;
while pointer+2*3+image_height*image_width*3-1<size(inputs,2)*size(inputs,1)
    temp = inputs(pointer+2*3:pointer+2*3+image_height*image_width*3-1);
    real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
    real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
    real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);
    
    pointer = pointer+2*3+image_height*image_width*3;
    i = i+1;
end


%inputs is 800x85000

length_of_time = ceil(triggers(end))+1;
upsampled_num_frames = length_of_time*120;

upsample_factor = round(refresh/(100/12));
bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
triggers = [triggers; triggers(end) + avg_bt_triggers];



for j = 1:size(spikes,2)
    
    
    frames_needed{j} = zeros(3,(length(triggers)-2)*100+120);
    a = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    frames_needed{j}(1,:) = a(1:(length(triggers)-2)*100+120);
    
    %     frames_needed{j}(1,:) = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    
    
    for i= 1:length(triggers)-1
        spacing = linspace(triggers(i), triggers(i+1),101);
        frames_needed{j}(2, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1);
    end
    
    for i = 1:size(frames_needed{j},2)-1
        frames_needed{j}(3,i) = sum(spikes{j} >= frames_needed{j}(2,i) & spikes{j} < frames_needed{j}(2,i+1));
    end
end
start_points = [1:10000:size(frames_needed{1},2) size(frames_needed{1},2)];
counter = 1;
sta_ref = datarun_ref.stas.stas{cell_indices_ref(1)};
[sig_stixels_ref] = significant_stixels(sta_ref, 'select', 'thresh', 'thresh', 3);


biggestBlob_stixs = ExtractNLargestBlobs(full(sig_stixels_ref), 1);
[matrix_subscript_i_eig, matrix_subscript_j_eig] = find(biggestBlob_stixs);
matrix_subscripts_eig = [matrix_subscript_i_eig, matrix_subscript_j_eig];

% compute the initial center of spatial fit
if size(matrix_subscripts_eig,1) > 1
    tmp = centroid(matrix_subscripts_eig);
    initial_center_point_x = tmp(2);
    initial_center_point_y = tmp(1);
else
    initial_center_point_x = matrix_subscripts_eig(2);
    initial_center_point_y = matrix_subscripts_eig(1);
end

initial_x = round(initial_center_point_x*stixel_size);
initial_y = round(initial_center_point_y*stixel_size);

        left = initial_x - 19;
        right = initial_x + 20;
        bottom = initial_y - 19;
        top = initial_y +20;

        if initial_x < 20
            left = 0;
        end
        if initial_x > 620
            right = 640;
        end
        if initial_y < 20;
            bottom = 0;
        end
        if initial_y > 300
            top = 320;
        end
        
% load('movie.mat')
movie = cell(1080,1);
for j = 1:length(start_points)-1
    for m = 1:10000/200
        fprintf(['%d out of ', num2str(length(start_points)-1), ' and %d out of ',num2str(10000/200), '\n' ], j, m)
        temp = load([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');


            movie{counter}.matrix = squeeze(temp.current_movie(left:right, bottom:top,2,:));% only use green channel


%         movie{counter}.matrix = squeeze(temp.current_movie(:,:,2,:));% only use green channel
        counter = counter+1;
        if counter > 1080
            break;
        end

    end

end

[spikesconcat, concat_fitmovie] = blocked_prep(blocked_spikes(1:1083), movie);


[STA, center] = STA_Test(spikesconcat, concat_fitmovie, 0);
[fittedGLM] = glm_fit(spikesconcat, concat_fitmovie, [20 20]);

plotfilters(fittedGLM)
