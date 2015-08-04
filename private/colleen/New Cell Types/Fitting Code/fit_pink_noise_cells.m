clear
close all

dataparam.date='2007-03-16-1';
dataparam.concatname='data003';

% Wrong Movie Information
% dataparam.file_name_wrong = [dataparam.date, '/', dataparam.concatname, '/data000', '/wrongMovie/wrongMovie'];
% dataparam.file_name_wrong = [dataparam.date, '/', dataparam.concatname, '/', 'wrongMovie/wrongMovie'];

% dataparam.mdf_file_wrong='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-4-0.48-22222.xml';

% Right Movie Information
dataparam.file_name_right = [dataparam.date, '/', 'data003-cr/', dataparam.concatname,'/', dataparam.concatname];
% dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname, '/','data000'];
%
dataparam.mdf_file_right='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-4-0.48-11111.xml';


% % Wrong Movie Information
% dataparam.mdf_file_wrong='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-22222.xml';
%
% % Right Movie Information
% dataparam.mdf_file_right='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
%
% %





% file_name = '2010-09-24-1/d05-36-norefit/data006-from-d05_36/fake_sta/data006-from-d05_36';
% mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-1-6-0.48-22222.xml';
%
% % Right Movie Information
% file_name2 = '2010-09-24-1/d05-36-norefit/data006-from-d05_cl36/real_sta/data006-from-d05_36';
% mdf_file2='/Volumes/Analysis/stimuli/white-noise-xml/BW-1-6-0.48-11111.xml';
%
%
% date='2010-09-24-1';
% concatname='data006-from-d05_36';



% fit parameters
fitparam.fit_n_one_filters = true;
fitparam.fit_n_two_filters = true;
fitparam.fit_surround_sd_scale = false;
fitparam.fit_surround = false;
fitparam.initial_n_one_filters =8;
fitparam.initial_n_two_filters = 8;

% fitparam.interpolate = false;
% fit the channels separately (useful for blue/green or SBCs)
fitparam.independent_fit = 0;
fitparam.num_frames = 50; % both have to be run with the name number of frames
% how many stixels per STA on average can be actually noise
fitparam.false_stixels =0.25;


dataparam.cell_type = {'unclassified'};





% list specific cell (1), or run for a whole cell type (0)
select_cells = 0;
% cell_specification = [218];%,861,1113,187,1307,1367,1669,1700,1787,1999,2133,2254,2315,2342,2417,2462]; %ON parasol
% dataparam.cell_specification = [17]; %OFF parasol

% cell_specification = [481];
%% END OF INPUT
dataparam.folder = dataparam.cell_type{1};
% file path to save data and pictures
dataparam.filepath=['/Users/colleen/Desktop/Fitting/',dataparam.date,'/',dataparam.concatname,'/'];
if ~exist([dataparam.filepath,dataparam.folder],'dir')
    mkdir([dataparam.filepath,dataparam.folder]);
end

% % Wrong Movie
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_wrong, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_wrong, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_wrong, '.sta'];

% Right Movie
datarun2.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun2.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun2.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];


%% Load Data1
% opt=struct('verbose',1,'load_params',0,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',0);
% opt.load_sta_params.save_rf = 1; % has to be set to one to load the data you need from vision
% opt.load_sta_params.frames =1:50;% have to input as a vector list of frames, not the number of frames total, counting backwards
% datarun=load_data(datarun,opt);

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
if select_cells == 0
    dataparam.cell_specification = datarun2.cell_types{cell_type_index}.cell_ids;
end



%% Movie for right STAs
% triggers=datarun2.triggers; %onsets of the stimulus presentation
%
% [mov,height,width,duration,refresh] = get_movie_ath(dataparam.mdf_file_right,...
%     triggers, 1,2);
%
% [mvi] = load_movie(dataparam.mdf_file_right, triggers);


% % check that stas have been loaded into datarun
% if ~isfield(datarun.stas, 'stas')
%     error('STAs are not contained in datarun, see LOAD_STA.M')
% end
%
% % initialize_output
% if ~isfield(datarun, 'matlab')
%     datarun = setfield(datarun, 'matlab', []);
% elseif ~isfield(datarun.matlab, 'sta_fits')
%     temp_cell = cell(length(datarun.cell_ids), 1);
%     datarun2.matlab.sta_fits = temp_cell;
% end


cell_indices = get_cell_indices(datarun2, dataparam.cell_specification);
num_rgcs = length(cell_indices);
variables = {'cell_specification', 'cell_indices', 'center_point_x', 'center_point_y', 'center_sd_x', 'center_sd_y', 'center_rotation_angle', 'color_weight_a', 'color_weight_b', 'color_weight_c', 'x_dim', 'y_dim', 'surround_sd_scale', 'surround_amp_scale', 'scale_one', 'scale_two', 'tau_one', 'tau_two','n_one_filters', 'n_two_filters','frame_number', 'rmse'};
parameters= zeros(num_rgcs,size(variables,2)); % want to save this

information{1} = dataset;
information{2} = variables;

paramPath = ['/Volumes/Analysis/', dataparam.date, '/', 'data003-cr/', dataparam.concatname,'/', dataparam.concatname, '.params'];

paramFile = edu.ucsc.neurobiology.vision.io.ParametersFile(paramPath);

for rgc = 1:num_rgcs
    
    % Vision STA for the wrong movie
    %     temp_sta = datarun.stas.stas{cell_indices(rgc)};
    
    % Determine the threshold that will result in the inputted false stixel
    % rate
    
    %     [threshold] = sig_stixels_threshold(temp_sta, fitparam.false_stixels);
    
    
    mark_params.select = 'thresh';
    mark_params.thresh = 2;%threshold;
    
    fprintf('fitting the STA for cell %d... \n', datarun2.cell_ids(cell_indices(rgc)))
    
    % Get the STA from the right movie
    temp_sta = datarun2.stas.stas{cell_indices(rgc)};
%     if fitparam.independent_fit
%         [temp_fit_params, sta, sta_temp, sig_stixels] = fit_sta_(temp_sta, 'fit_n_one_filters', fitparam.fit_n_one_filters,'fit_n_two_filters', fitparam.fit_n_two_filters, 'fit_surround_sd_scale', fitparam.fit_surround_sd_scale,'fit_surround', fitparam.fit_surround, 'initial_n_one_filters', fitparam.initial_n_one_filters,'initial_n_two_filters', fitparam.initial_n_two_filters, 'frame_number', fitparam.num_frames, 'mark_params', mark_params);
%         print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));
%         
%     else
        
        %         [temp_fit_params, sta, sig_stixels] = fit_sta(temp_sta, 'fit_n_one_filters', true, 'fit_n_two_filters', true, 'fit_surround_sd_scale', false, 'fit_surround', false, 'frame_number', num_frames, 'mark_params', mark_params);
        [temp_fit_params, sta,sig_stixels] = fit_sta_spatial(temp_sta, 'fit_n_one_filters', fitparam.fit_n_one_filters,'fit_n_two_filters', fitparam.fit_n_two_filters, 'fit_surround_sd_scale', fitparam.fit_surround_sd_scale,'fit_surround', fitparam.fit_surround, 'initial_n_one_filters', fitparam.initial_n_one_filters,'initial_n_two_filters', fitparam.initial_n_two_filters, 'frame_number', fitparam.num_frames, 'mark_params', mark_params,'biggest_blob', true);
        
        flip = strfind(lower(dataparam.cell_type{1}), 'off');
        if flip == 1
            plot_sta_fit(sta, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'off');
        else
            plot_sta_fit(sta, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'on');
            
        end
        
        suptitle({dataparam.dataset; sprintf('Fit for cell %d', datarun2.cell_ids(cell_indices(rgc)))})
        
        print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],['cell_',int2str(datarun2.cell_ids(cell_indices(rgc)))]));
        save_sig_stixels{rgc} = sig_stixels;
%     end
    
    % Plot result
    %  plot_sta_fit(sta_temp, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'off');
    
    
    % Add raw data to plot
    % hold on
    %         real_stix = significant_stixels(temp_sta);
    %       biggestBlob = ExtractNLargestBlobs(full(real_stix), 1);
    %     real_stix = biggestBlob;
    %     tc = time_course_from_sta(temp_sta, real_stix);
    %     norm_factor = max(abs(reshape(tc, 1, [])));
    %     tc = tc ./ norm_factor;
    %     hold on
    %     subplot(2,1,2)
    %         plot(tc(:,1), 'r', 'linewidth', 2)
    %         plot(tc(:,2), 'g', 'linewidth', 2)
    %         plot(tc(:,3), 'b', 'linewidth', 2)
    
    
    
    if isempty(temp_fit_params)
        temp_id = datarun.cell_ids(cell_indices(rgc));
        warn_message = ['cell ',num2str(temp_id), ' has no sig stixels and no fit'];
        warning(warn_message)
    end
    
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
    
    %
    %
    %
    % information{3} = parameters;
    % information{4} = datarun.matlab.sta_fits;
    % save([dataset,'-', cell_type{1}], 'information')
    
    tc = time_course_from_sta(sta, sig_stixels);
    norm_factor = max(abs(reshape(tc, 1, [])));
    tc = tc ./ norm_factor;
    %     if size(sta, 3) == 3
    %         plot(linspace(1,size(sta,4),size(tc,1)),tc(:,1), 'r')
    %         plot(linspace(1,size(sta,4),size(tc,1)),tc(:,2), 'g')
    %         plot(linspace(1,size(sta,4),size(tc,1)),tc(:,3), 'b')
    %     elseif size(sta_fit,3) == 1
    %         plot(linspace(1,size(sta,4),size(tc,1)), tc, 'k')
    %         hold off
    %     else
    %         error('dimensions of sta color is not recognized')
    %     end
    %     drawnow
%     stixel_size_vision = 29;
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaY', temp_fit_params.center_sd_x*5/64 );
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaX', temp_fit_params.center_sd_y*5/128 );
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'x0', temp_fit_params.center_point_x/128*5);
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', 5 - temp_fit_params.center_point_y/64*5);
    angle =tan(atan(-temp_fit_params.center_rotation_angle)*64/128); % transformation
    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'Theta', angle);
    if size(sta, 3) == 3
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', tc(:,1));
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,2));
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', tc(:,3));
    else
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,1));
  paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,1));
        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', tc(:,1));
%         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'blueness', 0); 
    end
    
    
    %         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaX', 0.2693);
    %         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaY', 0.9355);
    %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'x0', 3.0042);
    %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', 0.6847);
    %         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'Theta', 0.2738);
    %
    %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'nSpikes', 20000);
    
    
    
end
paramFile.close(1);

% return;

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


%% Look at the results
% figure
% suptitle({information{1}; [num2str(size(information{3},1)), cell_type{1}, ' cells']})
%
% %scale 1
% scale_one = information{3}(:,15);
% subplot(3,2,1)
% hist(scale_one);
% title('Scale One');
%
% scale_two = information{3}(:,16);
% subplot(3,2,2) ; hist(scale_two);
% title('Scale Two');
%
% tau_one = information{3}(:,17);
% subplot(3,2,3) ; hist(tau_one);
% title('Tau One');
%
% tau_two = information{3}(:,18);
% subplot(3,2,4) ; hist(tau_two);
% title('Tau Two');
%
% n = information{3}(:,19);
% subplot(3,2,5) ; hist(n);
% title('N Filters');
%
% area = information{3}(:,5).*information{3}(:,6)*pi;
% subplot(3,2,6) ; hist(area);
% title('RF Size');



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
    %     temp_rf = rf_from_sta(sta);
    %     imagesc(norm_image(temp_rf))
    %     hold on
    h(q) = plot_spatial_sd(all_params);
    set(h(q), 'DisplayName', 'Fitting Code');
    drawnow
end
axis equal

set(gca, 'xlim', [0, datarun2.stimulus.field_width]);
set(gca, 'ylim', [0, datarun2.stimulus.field_height]);
set(gca,'YDir','reverse');
%% compare mosaic to that from vision
% datarun=load_data('2010-09-24-0/data001-nwpca-cr');
% datarun=load_sta(datarun);
% datarun=load_params(datarun);
% datarun=load_neurons(datarun);
% datarun=set_polarities(datarun);
hold on
g = plot_rf_summaries(datarun2,{cell_type_index}, 'plot_fits',1,'label',0, 'foa', -1, 'clear', 0); %looking at the data type 8 is on large
children = get(g,'children');

legend_handle = legend([h(1), children(1)], {'Fitting Code', 'Vision'}, 'location', 'southeast');

title({['STA Fitting Code:  ', dataparam.cell_type{1}] ; dataset})
print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,dataparam.folder,'/'],[dataparam.cell_type{1}, ' Mosaic']));


%
% area  = parameters(:,5) .* parameters(:,6) .* pi;
% diameters = 2*[parameters(:,5) , parameters(:,6)];
% column1_larger = find(diameters(:,1)>diameters(:,2));
% column2_larger = find(diameters(:,1)<=diameters(:,2));
% max_diameters = [diameters(column1_larger,1); diameters(column2_larger,2)];
% disp('Mean RF diameter')
% mean(max_diameters)
% %% Calculate zero crossing
% zero_crossing = zeros(size(fit_tc,2),1);
% for i = 1:size(fit_tc,2)
%     cell_ = fit_tc{i}(:,2)*16.65; % units of ms
%
%     [x0,y0] = intersections(1:num_frames,cell_,1:num_frames, zeros(num_frames,1)); % find zero crossing
%     zero_crossing(i) = x0(1);
% end
% information{5} = zero_crossing;
% save([dataset,'-', cell_type{1}], 'information')
fitting_results = datarun2.matlab.sta_fits;
% [init_mean, init_sd, init_color_weight, init_scale, init_tau, init_n_filters] = computeAverageFittingParameters(fitting_results)
% save([sprintf('%s%s',[dataparam.filepath,dataparam.folder,'/']), 'results.mat'], 'fitting_results')


%
% cell_indices = get_cell_indices(datarun2, dataparam.cell_specification);
%
% paramPath = ['/Volumes/Analysis/', dataparam.date, '/', 'data003-cr/', dataparam.concatname,'/', dataparam.concatname, '.params'];
%
% paramFile = edu.ucsc.neurobiology.vision.io.ParametersFile(paramPath);
%
%
% for rgc= 1:1%length(datarun2.cell_ids)
%
%     stixel_size_vision = 29;
%    paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', 2.5);
%         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'x0', 3);
%
% %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', 5 -datarun2.vision.sta_fits{rgc}.mean(1));
% %         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'x0', datarun2.vision.sta_fits{rgc}.mean(2)*5);
%     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaX', datarun2.vision.sta_fits{rgc}.sd(1));
% %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'y0', (datarun2.stimulus.field_height -temp_fit_params.center_point_y)/stixel_size_vision );
%        paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'SigmaY', datarun2.vision.sta_fits{rgc}.sd(1));
%
%     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'Theta', datarun2.vision.sta_fits{rgc}.angle);
% %     if size(sta, 3) == 3
% %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', tc(:,1));
% %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,2));
% %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', tc(:,1));
% %     else
% %             paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'GreenTimeCourse', tc(:,1));
% % %                 paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'BlueTimeCourse', []);
% % %     paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'RedTimeCourse', []);
% %
% %                         paramFile.setCell(datarun2.cell_ids(cell_indices(rgc)), 'blueness', 0);
% %
% %
% %     end
%
%
% end
% paramFile.close(1);
