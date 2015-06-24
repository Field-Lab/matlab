% Find midget and parasol zero crossing and ratio of peaks

% load data


file_names  = {'2010-09-24-0/data001-nwpca/data001-nwpca', '2009-04-13-7/data000/data000', '2009-12-03-0/data003-nwpca/daat003/daat003', '2005-05-26-0/data001/data001', '2009-11-14-0/data000/data000', '2009-12-03-2/data000-nwpca/data000/data000', '2005-04-06-4/data000/data000', '2006-05-05-0/data000/data000', '2007-08-21-4/data000/data000', '2008-08-27-6/data000/data000', '2007-09-18-3/data000-mg/data000/data000', '2006-05-04-4/data000-nwpca/data000/data000', '2005-09-09-1/data000-nwpca/data000', '2008-03-25-4/data004/data004/data004' , '2005-05-26-7/data000/data000', '2005-08-03-0/data001/data001' };

% zero_crossing = zeros(size(file_names, 2),size(4,2));
% peak_ratio = zeros(size(file_names, 2),size(4,2));

    for files = 1:size(file_names, 2)
        file_name = file_names{files};
        datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
        datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
        datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
        
        %% Load Data1
        opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',0);
        opt.load_sta_params.save_rf = 1;
        opt.load_sta_params.frames =1:30;% have to input as a vector list of frames, not the number of frames total, counting backwards
        datarun=load_data(datarun,opt);
        
        cell_types = {'ON parasol', 'ON midget', 'OFF parasol', 'OFF midget'};
        
        
        for cell_type = 1:size(cell_types,2)
            if cell_type == 1 || cell_type == 2
                on = 1;
            else
                on = 0;
            end
            
            [indicies] = get_cell_indices(datarun, cell_types(cell_type));
            
            timecourses = zeros(length(indicies), 30);
            for i = 1:length(indicies)
                timecourses(i,:) = datarun.vision.timecourses(indicies(i)).g;
            end
            
            
            outlier_idx = abs(timecourses(:,27) - median(timecourses(:,27))) > 3*std(timecourses(:,27)); % Find outlier idx
            timecourses_good = timecourses(~outlier_idx,:);
            avg_timecourse = mean(timecourses_good,1)';
            
            
            
            time = [-29*datarun.stimulus.interval*8.33:datarun.stimulus.interval*8.33:0]';
            
            % figure
            % plot(time,avg_timecourse)
            % hold on
            % plot(time, zeros(size(time)), 'k')
            if on ~= 1
                avg_timecourse = -avg_timecourse;
            end
            
            [ind,t0] = crossing(avg_timecourse(13:30));
            t0= t0(1)+12;
            
            [Y,I] = max(avg_timecourse(13:30))
            peak_time = I(1)+ 12;
            
             [Y,I] = min(avg_timecourse(13:30))
            trough_time = I(1)+ 12;
            
            zero_crossing(files,cell_type) = (30-t0)*datarun.stimulus.interval*8.33;
            peak_ratio(files,cell_type) = abs(max(avg_timecourse))/abs(min(avg_timecourse));
            peak_to_zero(files, cell_type)  = peak_time - t0;
            trough_to_zero(files, cell_type)  = t0 - trough_time;
            peak_to_trough(files,cell_type) = peak_time - trough_time;
     

        end
        clear datarun
    end

