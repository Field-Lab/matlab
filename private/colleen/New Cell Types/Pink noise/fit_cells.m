clear
close all

dataparam.date='2015-08-17-1';
dataparam.concatname='d01-29-norefit/data026';


% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', 'data026'];
paramPath = ['/Volumes/Analysis/', dataparam.date, '/', dataparam.concatname, '/data026.params'];

dataparam.cell_type = {'unclassified'};

fitparam.num_frames = 50;
mark_params.thresh = 1.5; %threshold for significant stixels;
num_gauss = 0.5; % number of SD in fit


% list specific cell (1), or run for a whole cell type (0)
select_cells = 0;
if select_cells == 1
    dataparam.cell_specification = [31] %ON parasol
end

%% END OF INPUT
dataparam.folder = dataparam.cell_type{1};
% file path to save data and pictures
dataparam.filepath=['/Users/colleen/Desktop/Fitting/',dataparam.date,'/',dataparam.concatname,'/data026/'];
if ~exist([dataparam.filepath,dataparam.folder],'dir')
    mkdir([dataparam.filepath,dataparam.folder]);
end

% Right Movie
datarun2.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun2.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun2.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];



%% Load Data2
slashes = strfind(datarun2.names.rrs_neurons_path, '/');
dataset = datarun2.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun2=load_data(datarun2,opt);




%% Cell indicies

% Find the type of the inputted cells
cell_type_index= zeros(1,size(dataparam.cell_type,2));
for num_cell_types = 1:size(dataparam.cell_type,2)
    for i = 1:size(datarun2.cell_types,2)
        right_cell_type = strcmpi(datarun2.cell_types{i}.name, dataparam.cell_type{num_cell_types}); % case insensitive
        if right_cell_type == 1;
            cell_type_index(num_cell_types) = i;
            break
        end
        cell_type_index(num_cell_types) = 0;% couldn't find the right cell type
    end
    
end

% Set the cell_specification to all the cell of the inputted type
if select_cells ~= 1
    dataparam.cell_specification = datarun2.cell_types{cell_type_index}.cell_ids;
end


cell_indices = get_cell_indices(datarun2, dataparam.cell_specification);
num_rgcs = length(cell_indices);


paramFile = edu.ucsc.neurobiology.vision.io.ParametersFile(paramPath);

for rgc = 1:num_rgcs
    
    mark_params.select = 'thresh';
    
    fprintf('fitting the STA for cell %d... \n', datarun2.cell_ids(cell_indices(rgc)))
    
    % Get the STA from the right movie
    temp_sta = datarun2.stas.stas{cell_indices(rgc)};
    [final_fit_params, x,y] = fit_just_spatial(temp_sta, num_gauss, mark_params);
    h = final_fit_params(1);
    k = final_fit_params(2);
    a = num_gauss*final_fit_params(3);
    b = num_gauss*final_fit_params(4);
    angle = final_fit_params(5);
    A = final_fit_params(5);
    
    xhat = (x - h)*cos(angle) - (y-k)*sin(angle);
    yhat = (x - h)*sin(angle) + (y-k)*cos(angle);
    U = (xhat/a).^2 + (yhat/b).^2;
    F = A*exp(-U/2);
    
    fit_shaped = reshape(F, size(temp_sta,1),size(temp_sta,2));

    rf = rf_from_sta(temp_sta);
    if size(rf,3) > 1
        sta_one = rf(:,:,2);
    else
        sta_one =rf;
    end
    
    sta = temp_sta;

    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaY', b);
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaX', a);
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'x0', (h-0.5));
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', (size(sta,1) - (k))+0.5);

    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'Theta',  angle); % set all of them to 0 because of plotting problems
    sig_stixels = significant_stixels(double(sta), 'select', 'thresh', 'thresh',mark_params.thresh);
    if sum(full(sig_stixels)) == 0
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', nan(fitparam.num_frames,1));
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', nan(fitparam.num_frames,1));
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', nan(fitparam.num_frames,1));
        
    else
        
        tc = time_course_from_sta(sta,sig_stixels );
        norm_factor = max(abs(reshape(tc, 1, [])));
        tc = tc ./ norm_factor;
        
        if size(sta, 3) == 3
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', tc(:,1));
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,2));
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', tc(:,3));
        else
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', tc(:,1));
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,1));
            paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', tc(:,1));
            %         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'blueness', 0);
        end
    end
   
    
    
end
paramFile.close(1);
