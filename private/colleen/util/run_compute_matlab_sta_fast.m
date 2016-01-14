
clear
close all

dataparam.date='2008-11-12-1';
dataparam.concatname='data008-mg';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
dataparam.num_frames = 30;
dataparam.to_save = 1;
dataparam.scale = 10; %upsampling of STA
dataparam.save_location = ['/Users/colleen/Desktop/MATLAB_STAs/', dataparam.date, '/', dataparam.concatname, '/'];

% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/data008'];



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
opt.load_sta_params.frames = 1:dataparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);


% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
    vision_id = [1323] %ON parasol
else
    dataparam.cell_type = {'ON large 1'};
    cell_indices = get_cell_indices(datarun, dataparam.cell_type);
    vision_id=datarun.cell_ids(cell_indices);
end


[sta, timecourse, sig_stixels] = compute_matlab_sta_fast(datarun, dataparam, vision_id);
