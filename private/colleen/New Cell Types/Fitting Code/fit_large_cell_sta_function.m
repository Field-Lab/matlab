function fit_large_cell_sta_function(dataparam, fitparam)


%% END OF INPUT
dataparam.folder = dataparam.cell_type{1};
% file path to save data and pictures
if ~exist([dataparam.filepath,dataparam.folder],'dir')
    mkdir([dataparam.filepath,dataparam.folder]);
end

% Wrong Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_wrong, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_wrong, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_wrong, '.sta'];

% Right Movie
datarun2.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun2.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun2.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];


%% Load Data1
opt=struct('verbose',1,'load_params',0,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',0);
opt.load_sta_params.save_rf = 1; % has to be set to one to load the data you need from vision
opt.load_sta_params.frames =1:7;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

%% Load Data2
slashes = strfind(datarun2.names.rrs_neurons_path, '/');
dataset = datarun2.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun2=load_data(datarun2,opt);
if exist([dataparam.filepath, '/datarun.mat'])
	load([dataparam.filepath, '/datarun'], 'datarun2')
end



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
if cell_type_index == 0
    return;
end

% Set the cell_specification to all the cell of the inputted type
if dataparam.select_cells == 0
    dataparam.cell_specification = [];
    for i = 1:size(dataparam.cell_type,2)
        
        dataparam.cell_specification = [dataparam.cell_specification, datarun2.cell_types{cell_type_index(i)}.cell_ids];
    end
    
end



%% Movie for right STAs
triggers=datarun2.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(dataparam.mdf_file_right,...
    triggers, 1,2);

[mvi] = load_movie(dataparam.mdf_file_right, triggers);


% check that stas have been loaded into datarun
if ~isfield(datarun.stas, 'stas')
    error('STAs are not contained in datarun, see LOAD_STA.M')
end

% initialize_output
if ~isfield(datarun, 'matlab')
    datarun = setfield(datarun, 'matlab', []);
elseif ~isfield(datarun.matlab, 'sta_fits')
    temp_cell = cell(length(datarun.cell_ids), 1);
    datarun2.matlab.sta_fits = temp_cell;
end


cell_indices = get_cell_indices(datarun2, dataparam.cell_specification);
cell_ids=get_cell_ids(datarun2,dataparam.cell_specification);
num_rgcs = length(cell_indices);
variables = {'cell_specification', 'cell_indices', 'center_point_x', 'center_point_y', 'center_sd_x', 'center_sd_y', 'center_rotation_angle', 'color_weight_a', 'color_weight_b', 'color_weight_c', 'x_dim', 'y_dim', 'surround_sd_scale', 'surround_amp_scale', 'scale_one', 'scale_two', 'tau_one', 'tau_two','n_one_filters', 'n_two_filters','frame_number', 'rmse'};
parameters= zeros(num_rgcs,size(variables,2)); % want to save this

information{1} = dataset;
information{2} = variables;



for rgc = 1:num_rgcs
    
    % Vision STA for the wrong movie
    temp_sta = datarun.stas.stas{cell_indices(rgc)};
    
    % Determine the threshold that will result in the inputted false stixel
    % rate
    
    [threshold] = sig_stixels_threshold(temp_sta, fitparam.false_stixels);
    
    mark_params.select = 'thresh';
    mark_params.thresh = threshold;
    
    fprintf('fitting the STA for cell %d... \n', datarun2.cell_ids(cell_indices(rgc)))
    
    % Get the STA from the right movie
    temp_sta = datarun2.stas.stas{cell_indices(rgc)};
    
    % doesn't work well
    if fitparam.independent_fit
        [temp_fit_params, sta, sta_temp, sig_stixels] = fit_sta_(temp_sta, 'fit_n_one_filters', fitparam.fit_n_one_filters,'fit_n_two_filters', fitparam.fit_n_two_filters, 'fit_surround_sd_scale', fitparam.fit_surround_sd_scale,'fit_surround', fitparam.fit_surround, 'initial_n_one_filters', fitparam.initial_n_one_filters,'initial_n_two_filters', fitparam.initial_n_two_filters, 'fit_surround_amp_scale', fitparam.fit_surround_amp_scale, 'frame_number', fitparam.num_frames, 'mark_params', mark_params);
        print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));
        
    else
        [temp_fit_params, sta, sig_stixels] = fit_sta(temp_sta, 'fit_n_one_filters', fitparam.fit_n_one_filters,'fit_n_two_filters', fitparam.fit_n_two_filters, 'fit_surround_sd_scale', fitparam.fit_surround_sd_scale,'fit_surround', fitparam.fit_surround, 'initial_n_one_filters', fitparam.initial_n_one_filters,'initial_n_two_filters', fitparam.initial_n_two_filters, 'fit_surround_amp_scale', fitparam.fit_surround_amp_scale,'frame_number', fitparam.num_frames, 'mark_params', mark_params,'biggest_blob', false, 'verbose', false);
        
        flip = strfind(lower(dataparam.cell_type{1}), 'off');
        if flip == 1
            plot_sta_fit(sta, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'off');
        else
            plot_sta_fit(sta, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'on');
            
        end
        
        suptitle({dataparam.dataset; sprintf('Fit for cell %d', datarun2.cell_ids(cell_indices(rgc)))})
        
        print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],['cell_',int2str(datarun2.cell_ids(cell_indices(rgc)))]));
        save_sig_stixels{rgc} = sig_stixels;
    end
    
 
    
    
    if isempty(temp_fit_params)
        temp_id = datarun.cell_ids(cell_indices(rgc));
        warn_message = ['cell ',num2str(temp_id), ' has no sig stixels and no fit'];
        warning(warn_message)
    end
    
%     if temp_fit_params.center_rotation_angle < 0 
%         temp_fit_params.center_rotation_angle = 2*pi - temp_fit_params.center_rotation_angle;
%     end
    
    datarun2.matlab.sta_fits{cell_indices(rgc)} = temp_fit_params;
    datarun2.matlab.sta_fits{cell_indices(rgc)}.sig_stixels = sig_stixels;
    
    
    
    parameters(rgc,1) = dataparam.cell_specification(rgc);
    parameters(rgc,2) = cell_indices(rgc);
    parameters(rgc,3) = datarun2.matlab.sta_fits{cell_indices(rgc)}.center_point_x;
    parameters(rgc,4) = datarun2.matlab.sta_fits{cell_indices(rgc)}.center_point_y;
    parameters(rgc,5) = datarun2.matlab.sta_fits{cell_indices(rgc)}.center_sd_x;
    parameters(rgc,6) = datarun2.matlab.sta_fits{cell_indices(rgc)}.center_sd_y;
    parameters(rgc,7) = datarun2.matlab.sta_fits{cell_indices(rgc)}.center_rotation_angle;
    parameters(rgc,8) = datarun2.matlab.sta_fits{cell_indices(rgc)}.color_weight_a;
    parameters(rgc,9) = datarun2.matlab.sta_fits{cell_indices(rgc)}.color_weight_b;
    parameters(rgc,10) = datarun2.matlab.sta_fits{cell_indices(rgc)}.color_weight_c;
    parameters(rgc,11) = datarun2.matlab.sta_fits{cell_indices(rgc)}.x_dim;
    parameters(rgc,12) = datarun2.matlab.sta_fits{cell_indices(rgc)}.y_dim;
    parameters(rgc,13) = datarun2.matlab.sta_fits{cell_indices(rgc)}.surround_sd_scale;
    parameters(rgc,14) = datarun2.matlab.sta_fits{cell_indices(rgc)}.surround_amp_scale;
    parameters(rgc,15) = datarun2.matlab.sta_fits{cell_indices(rgc)}.scale_one;
    parameters(rgc,16) = datarun2.matlab.sta_fits{cell_indices(rgc)}.scale_two;
    parameters(rgc,17) = datarun2.matlab.sta_fits{cell_indices(rgc)}.tau_one;
    parameters(rgc,18) = datarun2.matlab.sta_fits{cell_indices(rgc)}.tau_two;
    parameters(rgc,19) = datarun2.matlab.sta_fits{cell_indices(rgc)}.n_one_filters;
    parameters(rgc,20) = datarun2.matlab.sta_fits{cell_indices(rgc)}.n_two_filters;
    
    parameters(rgc,21) = datarun2.matlab.sta_fits{cell_indices(rgc)}.frame_number;
    parameters(rgc,22) = datarun2.matlab.sta_fits{cell_indices(rgc)}.rmse;
    
    sta = datarun2.stas.stas{cell_indices(rgc)};
    parameters_per_cell = datarun2.matlab.sta_fits{cell_indices(rgc)};
    save(sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],['cell_',int2str(datarun2.cell_ids(cell_indices(rgc)))]), 'dataparam', 'fitparam', 'parameters_per_cell', 'sta' )
    save([dataparam.filepath, '/datarun'], 'datarun2')

    
end




figure
for q = 1:num_rgcs
    temp_sta = datarun2.stas.stas{cell_indices(q)};
    temp_fit_params = datarun2.matlab.sta_fits{cell_indices(q)};
    fit_tc{q} =  plot_fit_timecourses(temp_sta, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, temp_fit_params.fit_params, temp_fit_params.fixed_params, save_sig_stixels{rgc}, 1); % plot_raw = 1
    hold on
    
end
title({['Fits of ' num2str(length(dataparam.cell_specification)), ' ', dataparam.cell_type{1}, ' Cells']; dataset})
xlabel('Time')
ylabel('Amplitude')

print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],[dataparam.cell_type{1}, ' Timecourses']));

%% Plot mosaic
figure
for q = 1:num_rgcs
    temp_sta = datarun2.stas.stas{cell_indices(q)};
    temp_fit_params = datarun2.matlab.sta_fits{cell_indices(q)};
    fit_indices = temp_fit_params.fit_indices;
    fixed_indices = temp_fit_params.fixed_indices;
    fit_params = temp_fit_params.fit_params;
    fixed_params = temp_fit_params.fixed_params;
    
    all_params(fit_indices) = fit_params;
    all_params(fixed_indices) = fixed_params;
    
    % get sta fit
    sta_fit = sta_fit_function(all_params);
    % spatial fit
    
    
    hold on
    h(q) = plot_spatial_sd(all_params);
    set(h(q), 'DisplayName', 'Fitting Code');
    drawnow
end
axis equal

set(gca, 'xlim', [0, datarun2.stimulus.field_width]);
set(gca, 'ylim', [0, datarun2.stimulus.field_height]);
set(gca,'YDir','reverse');
%% compare mosaic to that from vision
% 
% hold on
% g = plot_rf_summaries(datarun2,{cell_type_index}, 'plot_fits',1,'label',0, 'foa', -1, 'clear', 0); %looking at the data type 8 is on large
% children = get(g,'children');
% 
% legend_handle = legend([h(1), children(1)], {'Fitting Code', 'Vision'}, 'location', 'southeast');
% 
% title({['STA Fitting Code:  ', dataparam.cell_type{1}] ; dataset})
% print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],[dataparam.cell_type{1}, ' Mosaic']));

fitting_results = datarun2.matlab.sta_fits;

%% Write to Vision File

    
paramFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun2.names.rrs_params_path);
for i = 1:length(cell_ids)
paramFile.setCell(cell_ids(i), 'x0', datarun2.matlab.sta_fits{cell_indices(i)}.center_point_x-0.5)
paramFile.setCell(cell_ids(i), 'y0', datarun2.stimulus.field_height - datarun2.matlab.sta_fits{cell_indices(i)}.center_point_y+0.5)

paramFile.setCell(cell_ids(i), 'SigmaX', datarun2.matlab.sta_fits{cell_indices(i)}.center_sd_y)

paramFile.setCell(cell_ids(i), 'SigmaY', datarun2.matlab.sta_fits{cell_indices(i)}.center_sd_x)
paramFile.setCell(cell_ids(i), 'Theta', (2*pi - datarun2.matlab.sta_fits{cell_indices(i)}.center_rotation_angle))

end

paramFile.close(1);
% ind = get_cell_indices(datarun2, 'ON large 1');
% for i = 1:length(ind)
%     large_dia_on(i) = geomean([datarun2.matlab.sta_fits{ind(i)}.center_sd_x,datarun2.matlab.sta_fits{ind(i)}.center_sd_y])*2./(datarun2.matlab.sta_fits{ind(i)}.x_dim/640).*5.5;
% end
% 
% ind = get_cell_indices(datarun2, 'OFF large 1');
% for i = 1:length(ind)
%     large_dia_off(i) = geomean([datarun2.matlab.sta_fits{ind(i)}.center_sd_x,datarun2.matlab.sta_fits{ind(i)}.center_sd_y])*2./(datarun2.matlab.sta_fits{ind(i)}.x_dim/640).*5.5;
% end
