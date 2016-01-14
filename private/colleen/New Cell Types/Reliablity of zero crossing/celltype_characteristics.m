% cell properties


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_opts.date='2007-05-01-2'; % one slash at the end
run_opts.concatname='data000'; % Name (or modified name) of run, no slashes

% Sometimes the data has two versions of the concate name
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data000'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
%
% Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';


run_opts.save_location_root = '/Users/colleen/Desktop/Cell Characteristics/';
% Where to save
run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, ];


% Cell specification can be one cell type or multiple in a cell array.
% Use the same spelling/capitalization as the vision params file
% OFF large-2 in vision = OFF large 2
run_opts.cell_specification = {'ON large 2'};

if isempty(strfind(run_opts.cell_specification{1}, 'large'))
    run_opts.fitted = 0;
else
    run_opts.fitted = 0;
end


run_opts.false_stixels = 0.5; % increase if stixels are size 32
run_opts.num_frames = 30;
run_opts.frames_past_zero = 0; %on the vision interface

if run_opts.fitted ==1;
    load_path = '/Users/colleen/Desktop/Fitting/';
    full_path = [load_path, run_opts.date, '/', run_opts.concatname, '/', run_opts.cell_specification{1}];
    s=dir(fullfile(full_path,'*.mat'));
    
    for i = 1:size(s,1)
        fitted_cells(i) =  load([full_path, '/', s(i).name]);
        
    end
    for i = 1:size(s,1)
        period_location = strfind(s(i).name, '.');
        fitted_cells(i).name = str2num(s(i).name(6:period_location-1)); % add in the vision id of the cell
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
datarun=load_data(datarun,opt);

refresh = datarun.stimulus.interval/120*1000;


if run_opts.fitted == 1
    for i =1:size(s,1)
        [indicies(i)] = get_cell_indices(datarun, fitted_cells(i).name);
    end
    fig = figure;
    for i = 1:length(indicies)
        figure(fig);
        sta = fitted_cells(i).sta;
        all_params(fitted_cells(i).parameters_per_cell.fit_indices) = fitted_cells(i).parameters_per_cell.fit_params;
        all_params(fitted_cells(i).parameters_per_cell.fixed_indices) = fitted_cells(i).parameters_per_cell.fixed_params;
        
        % get sta fit
        sta_fit = sta_fit_function(all_params);
        
        temp_stix = fitted_cells(i).parameters_per_cell.sig_stixels;
        % temp_stix = significant_stixels(sta); %changed from 3.5 to 3.0
        % biggestBlob = ExtractNLargestBlobs(full(temp_stix), 1);
        % temp_stix = biggestBlob;
        % temp_stix = sig_stixels
        % fit_tc = time_course_from_sta(sta_fit, temp_stix);
        fit_tc_full = time_course_from_sta(sta_fit, temp_stix);
        
        % norm_factor = max(abs(reshape(fit_tc_full, 1, [])));
        % fit_tc{i} = fit_tc_full ./ norm_factor;
        fit_tc{i} = fit_tc_full;
        
        if size(fit_tc{i},2) > 1
            timecourse_green(i,:) = fit_tc{i}(:,2);
        else
            timecourse_green(i,:) = fit_tc{i}(:,1);
        end
        
        
        if size(sta_fit, 3) == 3
            plot(linspace(1,size(sta,4),size(fit_tc{i},1)), fit_tc{i}(:,1), 'r')
            hold on
            plot(linspace(1,size(sta,4),size(fit_tc{i},1)),fit_tc{i}(:,2), 'g')
            t_line = plot(linspace(1,size(sta,4),size(fit_tc{i},1)),fit_tc{i}(:,3), 'b');
        elseif size(sta_fit, 3) == 1
            t_line = plot(linspace(1,size(sta,4),size(fit_tc{i},1)), fit_tc{i}, 'k');
            hold on
        else
            error('dimensions of sta color is not recognized')
        end
        
        
        if abs(min(timecourse_green(i,:))) > abs(max(timecourse_green(i,:)))
            flip = 1;
        else
            flip =0;
        end
        
        
        % get sta fit
        % spatial fit
        figure
        
        temp_rf = rf_from_sta(sta, 'sig_stixels', temp_stix);
        if flip == 1
            sta_img = imagesc(norm_image(-temp_rf));
            hold on
            c = plot_spatial_sd(all_params);
            set(c, 'color', 'k')
        else
            sta_img = imagesc(norm_image(temp_rf));
            hold on
            c = plot_spatial_sd(all_params);
            set(c, 'color', 'k')
        end
        
        
        axis image
        drawnow
        
        sta_axes{i} = get(sta_img, 'Parent');
        
    end
else
    indicies = get_cell_indices(datarun, run_opts.cell_specification{1});
    %     indicies=get_cell_ids(datarun,run_opts.cell_specification{1});
    fig = figure;
    for i = 1:length(indicies)
        index = indicies(i);
        figure(fig);
        
        sta = datarun.stas.stas{index};
        %         all_params(fitted_cells(i).parameters_per_cell.fit_indices) = fitted_cells(i).parameters_per_cell.fit_params;
        %         all_params(fitted_cells(i).parameters_per_cell.fixed_indices) = fitted_cells(i).parameters_per_cell.fixed_params;
        
        % get sta fit
        %         sta_fit = sta_fit_function(all_params);
        
        %         temp_stix = fitted_cells(i).parameters_per_cell.sig_stixels;
        % temp_stix = significant_stixels(sta); %changed from 3.5 to 3.0
        % biggestBlob = ExtractNLargestBlobs(full(temp_stix), 1);
        % temp_stix = biggestBlob;
        % temp_stix = sig_stixels
        % fit_tc = time_course_from_sta(sta_fit, temp_stix);
        try
            [threshold] = sig_stixels_threshold(sta, run_opts.false_stixels);
        catch
            threshold = 2;
        end
        
        [sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', threshold);
        fit_tc_full = time_course_from_sta(sta, sig_stixels);
        
        % norm_factor = max(abs(reshape(fit_tc_full, 1, [])));
        % fit_tc{i} = fit_tc_full ./ norm_factor;
        fit_tc{i} = fit_tc_full;
        
        if ~isempty(fit_tc{i})
            if size(fit_tc{i},2) > 1
                timecourse_green(i,:) = fit_tc{i}(:,2);
            else
                timecourse_green(i,:) = fit_tc{i}(:,1);
            end
        

            if size(sta, 3) == 3
                plot(linspace(1,size(sta,4),size(fit_tc{i},1)), fit_tc{i}(:,1), 'r')
                hold on
                plot(linspace(1,size(sta,4),size(fit_tc{i},1)),fit_tc{i}(:,2), 'g')
                t_line = plot(linspace(1,size(sta,4),size(fit_tc{i},1)),fit_tc{i}(:,3), 'b');
            elseif size(sta, 3) == 1
                t_line = plot(linspace(1,size(sta,4),size(fit_tc{i},1)), fit_tc{i}, 'k');
                hold on
            else
                error('dimensions of sta color is not recognized')
            end
                
        
            if abs(min(timecourse_green(i,:))) > abs(max(timecourse_green(i,:)))
                flip = 1;
            else
                flip =0;
            end
end
        
        figure
        temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
        
                if isempty(fit_tc{i})
flip = 0;
                end
                
                    
        if flip == 1
            sta_img = imagesc(norm_image(-temp_rf));
            hold on
            
            
        else
            sta_img = imagesc(norm_image(temp_rf));
            hold on
        end
        the_fit = datarun.stas.fits{index};
        ctr = the_fit.mean;
        rad = the_fit.sd;
        [X,Y] = drawEllipse([ctr rad the_fit.angle]);
        plot(X,Y,'Color','k');
        hold on
        set(gca, 'ydir', 'reverse')
        
        axis image
        drawnow
        
        sta_axes{i} = get(sta_img, 'Parent');
        
    end
end

t_axes = get(t_line, 'Parent');

outlier_idx = abs(timecourse_green(:,27) - median(timecourse_green(:,27))) > 3*std(timecourse_green(:,27)); % Find outlier idx

timecourses_good = timecourse_green(~outlier_idx,:);
avg_timecourse = mean(timecourses_good,1)';
if abs(min(avg_timecourse)) > abs(max(avg_timecourse))
    avg_timecourse_orient  = -avg_timecourse;
else
    avg_timecourse_orient = avg_timecourse;
end

time = -(run_opts.num_frames-run_opts.frames_past_zero-1)*refresh:refresh:refresh + run_opts.frames_past_zero - 1;

[~,t0] = crossing(avg_timecourse_orient);
[~,ind] = min(abs(t0 - run_opts.num_frames/2)); % find the frame that closest to the middle of num_frames
t0 = t0(ind);
zc_time = abs(t0*refresh + time(1)- refresh);

[peak , peak_ind] = max(avg_timecourse_orient);
peak_time = abs(peak_ind*refresh + time(1)- refresh);
amplitude = peak;

[trough , trough_ind] = min(avg_timecourse_orient);
trough_time = abs(trough_ind*refresh + time(1)- refresh);

% num_spikes = 0;
if run_opts.fitted == 1
    gmean = nan(length(indicies),1);
    for i = 1:length(indicies)
        num_spikes(i) = length(datarun.spikes{indicies(i)});
        gmean(i) = geomean([fitted_cells(i).parameters_per_cell.center_sd_x*datarun.stimulus.stixel_width*5.5,fitted_cells(i).parameters_per_cell.center_sd_y*datarun.stimulus.stixel_width*5.5]); % 5.5 um/pixel
        
        
    end
else
    gmean = nan(length(indicies),1);
    for i = 1:length(indicies)
        num_spikes(i) = length(datarun.spikes{indicies(i)});
        gmean(i) = geomean([datarun.stas.fits{indicies(i)}.sd(1)*datarun.stimulus.stixel_width*5.5,datarun.stas.fits{indicies(i)}.sd(2)*datarun.stimulus.stixel_width*5.5]); % 5.5 um/pixel
        
        
    end
end

RF_diameter = gmean*2; % want the diameter instead of the radius
FR = num_spikes/datarun.duration;

% Plot mosaic from fits
if run_opts.fitted == 1
    figure
    for i = 1:length(indicies)
        temp_sta = fitted_cells(i).sta;
        fit_indices = fitted_cells(i).parameters_per_cell.fit_indices;
        fixed_indices = fitted_cells(i).parameters_per_cell.fixed_indices;
        fit_params = fitted_cells(i).parameters_per_cell.fit_params;
        fixed_params = fitted_cells(i).parameters_per_cell.fixed_params;
        
        all_params(fit_indices) = fit_params;
        all_params(fixed_indices) = fixed_params;
        
        % get sta fit
        sta_fit = sta_fit_function(all_params);
        % spatial fit
        
        
        hold on
        
        h(i) = plot_spatial_sd(all_params);
        set(h(i), 'color', 'k');
        drawnow
    end
    
    axis equal
    
    set(gca, 'xlim', [0, datarun.stimulus.field_width]);
    set(gca, 'ylim', [0, datarun.stimulus.field_height]);
    set(gca,'YDir','reverse');
    mosaic_axes = get(h(i), 'Parent');
    
    
else
    
    h = plot_rf_summaries(datarun, run_opts.cell_specification);
    
    
    axis equal
    
    set(gca, 'xlim', [0, datarun.stimulus.field_width]);
    set(gca, 'ylim', [0, datarun.stimulus.field_height]);
    set(gca,'YDir','reverse');
    mosaic_axes = h;
%     mosaic_axes = get(h, 'Parent');
end

% characterization structure;

prop.run_opts = run_opts;
prop.refresh = refresh;
prop.avg_timecourse = avg_timecourse;
prop.time = time;
prop.zc_time = zc_time;
prop.peak_time = peak_time;
prop.trough_time = trough_time;
prop.amplitude = amplitude;
prop.RF_diameter = RF_diameter;
prop.FR = FR;
prop.mosaic_axes = mosaic_axes;
prop.t_axes = t_axes;
prop.sta_axes = sta_axes;

% Make the directory to save in if it doesn't exist
if ~exist(run_opts.filepath)
    mkdir(run_opts.filepath)
end


save([run_opts.filepath, run_opts.cell_specification{1}], 'prop')