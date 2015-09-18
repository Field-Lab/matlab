% compute axon velocity for various cell types
clear
close all

dataparam.date='2007-03-02-0';
dataparam.concatname='d00_01_03_11_14_15';
% dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-6-0.48-11111.xml';
dataparam.num_frames = 30;
dataparam.to_save = 1;
% dataparam.scale = 10; %upsampling of STA
dataparam.save_location = ['/Users/colleen/Desktop/EIs/', dataparam.date, '/', dataparam.concatname, '/'];

% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/data003-from-data000_data001_data003_data011_data014_data015/data003-from-data000_data001_data003_data011_data014_data015'];

if ~exist(dataparam.save_location)
    mkdir(dataparam.save_location);
end


%% END OF INPUT
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];
datarun.names.rrs_ei_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.ei'];



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

dataparam.cell_types = {'ON parasol nc0', 'ON midget nc1', 'ON large 1', 'ON large 2'};


for cell_type = 1:length(dataparam.cell_types)
    
    select_cells = 0;
    if select_cells == 1
        vision_id = [1323]; %ON parasol
    else
        cell_indices = get_cell_indices(datarun, dataparam.cell_types{cell_type});
        vision_id=datarun.cell_ids(cell_indices);
    end

    datarun = load_ei(datarun, vision_id, 'array_type', 512);

    
    for i  = 1:length(vision_id)

        [v{cell_type}(i),c]=ei_axon_speed(datarun,vision_id(i), 'threshold', 3);
    end

end

%% Plot after it has run several times
figure; 
cmap = hsv(size(v,2));
counter = 1;
bins = 0:0.2:2;
for i  = 1:size(v,2)
    bins = 0.5:0.2+0.01*i:2;
    if sum(isnan(v{i})) ~= length(v{i})
        [n1, x_v] = hist(v{i},bins);
        hold on
        bar(x_v,n1, 'FaceColor', cmap(i,:)); 
        legend_entries{counter} = dataparam.cell_types{i};
        counter = counter +  1;
    else
        legend_entries{i} = '';
    end
    
end
axis tight
legend(legend_entries, 'location', 'northwest')
title({'Axon Velocity';dataparam.date;'data003'})
xlabel('Axon Velocity');
ylabel('Number of cells');

hgexport(gcf, [dataparam.save_location, 'data003']) 