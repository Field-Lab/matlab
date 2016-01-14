clear
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.sta';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
indicies = get_cell_indices(datarun, 544);
sta = datarun.stas.stas{indicies};
 

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 2);
figure
temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
image = norm_image(temp_rf);


fit_tc = time_course_from_sta(sta, sig_stixels);
        
        % norm_factor = max(abs(reshape(fit_tc_full, 1, [])));
        % fit_tc{i} = fit_tc_full ./ norm_factor;
        

            if size(sta, 3) == 3
                plot(linspace(1,size(sta,4),size(fit_tc,1)), fit_tc(:,1), 'r')
                hold on
                plot(linspace(1,size(sta,4),size(fit_tc,1)),fit_tc(:,2), 'g')
                t_line = plot(linspace(1,size(sta,4),size(fit_tc,1)),fit_tc(:,3), 'b');
            elseif size(sta, 3) == 1
                t_line = plot(linspace(1,size(sta,4),size(fit_tc,1)), fit_tc, 'k');
                hold on
            else
                error('dimensions of sta color is not recognized')
            end
        
                