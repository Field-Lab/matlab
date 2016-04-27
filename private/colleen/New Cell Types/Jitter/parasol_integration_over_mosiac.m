% parasol integration area over the mosaic 

clear
close all
% dbstop if error
dataparam.date='2016-02-17-6';
dataparam.concatname='data023';
dataparam.jitter_concatname='data026';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-33333-119.5.xml';
dataparam.stixel_size = 8;
dataparam.seed = 33333;
fitparam.num_frames = 20;
frame_width = 640/dataparam.stixel_size;
frame_height = 320/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;

num_colors =3;
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];

jitter_id = 6513;
load(['/Volumes/Lab/Users/crhoades/Jitter/', dataparam.date, '/', dataparam.jitter_concatname, '/Cell ', num2str(jitter_id), '.mat'])


% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
   dataparam.cell_specification = ['ON parasol']; %ON parasol

end
dataparam.cell_type = {'all'};
%% END OF INPUT


% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];


%% Load Data2

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_indices = get_cell_indices(datarun, dataparam.cell_specification);
total_sta = zeros(size(datarun.stas.stas{1}));
for i =1:length(cell_indices)
    total_sta = total_sta + double(datarun.stas.stas{cell_indices(i)});
    
end

