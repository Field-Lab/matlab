% runJitterAnalysis
% git hash 75ce3d76b7b20220cd94d83eb96e8e3ec90dd20e

clear
close all
% dbstop if error
dataparam.date='2016-02-17-6/data026_cf_split/edited';
dataparam.concatname='data026_cf_split';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-2-0.48-22222-119.5.xml';
dataparam.stixel_size = 16;
dataparam.seed = 22222;
fitparam.num_frames = 30;
frame_width = 640/dataparam.stixel_size;
frame_height = 320/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;

num_colors =3;
% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];

% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
%        dataparam.cell_specification = [184] %ON parasol
%   dataparam.cell_specification = [241 440 481 485 545 813 816 1386 1441 1534 1759 1786 1788 1936 1937 1940 2014 2103 2137 2162 2163 2165 2166 2167 2331 2405 2431 3048 3288 3436 3816 3888 3903 3907 3908 4113 4280 4937 4981 5013 5071 5522 5915 6110 6211 6336 6337 6512 6513 6529 7054 7460 7461 7462 ] %ON parasol
  dataparam.cell_specification = [1021 1022 1023 1027 1028 1741 1742 1743 1744 2751 2752 2753 3066 3067 3068 3466 3468 3469 3470 4296 4297 4298 4299 4300 4327 4328 4329 4330 4332 5026 5027 5028 5176 5177 5178 5179 5180 5536 5537 5538 5539] %ON parasol
 %  dataparam.cell_specification = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol
 %  dataparam.cell_specification = [185 1217 1221 1938 1942 3037 3904 3908 4893 4982 5523 5526 6094 6533 7056 7475 481 5061 3618 3811 5916 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 394 1218 1744 1970 2420 3302 4053 4097 4473 5987 6063 6617] %ON parasol
%dataparam.cell_specification = [481 812 1535 2102 2162 2987 3289 3812 4322 4939 6332 6512 1938 2012 3169 3903 4052 4981 5012 6541 7054 7458 6212 436 544 2432 3319 5074 1786 4549];
   %dataparam.cell_specification = [184 ];
    %dataparam.cell_specification = [481 814 1477 1986 2167 2168 2991 3289 4326 4941 6513 6517 1441 1786 2431 3436 5026 6212 545 6484 185 1217 1221 1938 1942 2012 3037 3904 3908 4112 4893 4982 5523 5526 6094 6533 7056 7475 5061 3618 3811] %ON parasol
 dataparam.cell_specification = [1021 ];   
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

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
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

    jitterX = load('jitterX_16.mat');
    jitterX = jitterX.jitterX;
        jitterY = load('jitterY_16.mat');
    jitterY = jitterY.jitterY;
for j = 1:length(cell_indices)
    
    
    spikes{j}=  datarun.spikes{cell_indices(j)};

end

%  state = Init_RNG_JavaStyle(seed);
%     jitterX = nan(duration,1);
%     jitterY = nan(duration,1);
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
%     save ('jitterX_8', 'jitterX');
%     save ('jitterY_8', 'jitterY');
    



[sta] = compute_jitter_sta_20160217(datarun, dataparam.mdf_file, fitparam.num_frames, spikes, jitterX, jitterY, stixel_size, num_colors, dataparam);
for i = 1:size(spikes,2)
    temp = sta{i};
    if ~exist(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/', dataparam.concatname])
        mkdir(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/',dataparam.concatname]);
    end
    save(['/Volumes/Lab/Users/crhoades/Jitter/',dataparam.date,'/', dataparam.concatname, '/Cell ', num2str(cell_ids(i))], 'temp')
end

toc

[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);

