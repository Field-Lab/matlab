% Find midget and parasol zero crossing and ratio of peaks

% load data


% file_names  = {'2010-09-24-0/data001-nwpca/data001-nwpca', '2009-04-13-7/data000/data000', '2009-12-03-0/data003-nwpca/daat003/daat003', '2005-05-26-0/data001/data001', '2009-11-14-0/data000/data000', '2009-12-03-2/data000-nwpca/data000/data000', '2005-04-06-4/data000/data000', '2006-05-05-0/data000/data000', '2007-08-21-4/data000/data000', '2008-08-27-6/data000/data000', '2007-09-18-3/data000-mg/data000/data000', '2006-05-04-4/data000-nwpca/data000/data000', '2005-09-09-1/data000-nwpca/data000', '2008-03-25-4/data004/data004/data004' , '2005-05-26-7/data000/data000', '2005-08-03-0/data001/data001' };
file_names  = {'2009-12-03-0/data003-nwpca/daat003/daat003', '2005-05-26-0/data001/data001', '2009-12-03-2/data000-nwpca/data000/data000', '2005-04-06-4/data000/data000', '2006-05-05-0/data000/data000', '2007-08-21-4/data000/data000', '2008-08-27-6/data000/data000', '2007-09-18-3/data000-mg/data000/data000', '2006-05-04-4/data000-nwpca/data000/data000', '2005-09-09-1/data000-nwpca/data000', '2008-03-25-4/data004/data004/data004' , '2005-05-26-7/data000/data000', '2005-08-03-0/data001/data001' };

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
            
            [ind,t0] = crossing(avg_timecourse(7:27));
            t0= t0(end)+8;
            
            [Y,I] = max(avg_timecourse(7:27))
            peak_time = I(end)+ 8;
            
             [Y,I] = min(avg_timecourse(7:27))
            trough_time = I(end)+ 8;
           
            zero_crossing(files,cell_type) = (30-t0)*datarun.stimulus.interval*8.33;
            peak_ratio(files,cell_type) = abs(max(avg_timecourse))/abs(min(avg_timecourse));
            peak_to_zero(files, cell_type)  = zero_crossing(files,cell_type) - (30-peak_time)*datarun.stimulus.interval*8.33;
            trough_to_zero(files, cell_type)  =(30-trough_time)*datarun.stimulus.interval*8.33 - zero_crossing(files,cell_type);
            peak_to_trough(files,cell_type) = (30-trough_time)*datarun.stimulus.interval*8.33 - (30-peak_time)*datarun.stimulus.interval*8.33;
            biphasic_index(files, cell_type) = abs(max(avg_timecourse)/min(avg_timecourse));

        end
        clear datarun
    end

zero_cross_ON = zero_crossing(:,2)./zero_crossing(:,1);    
figure; hist(zero_cross_ON)
xlabel('ON Parasol Zero Crossing/ON Midget Zero Crossing')
ylabel('Number of datasets')
title('Zero Crossing ON Cells')




zero_cross_OFF = zero_crossing(:,4)./zero_crossing(:,3)    
figure; hist(zero_cross_OFF)
xlabel('OFF Parasol Zero Crossing/OFF Midget Zero Crossing')
ylabel('Number of datasets')
title('Zero Crossing OFF Cells')

biphasic_index_ON = biphasic_index(:,2)./biphasic_index(:,1);    
figure; hist(biphasic_index_ON)
xlabel('ON Parasol Biphasic Index/ON Midget Biphasic Index')
ylabel('Number of datasets')
title('Biphasic Index ON Cells')

biphasic_index_OFF = biphasic_index(:,4)./biphasic_index(:,3);    
figure; hist(biphasic_index_OFF)
xlabel('OFF Parasol Biphasic Index/OFF Midget Biphasic Index')
ylabel('Number of datasets')
title('Biphasic Index OFF Cells')

figure
plot(ones(size(biphasic_index)), biphasic_index(:,1), 'ob')
hold on
plot(2*ones(size(biphasic_index)), biphasic_index(:,2), 'ob')
plot(3*ones(size(biphasic_index)), biphasic_index(:,3), 'ob')
plot(4*ones(size(biphasic_index)), biphasic_index(:,4), 'ob')

xlabel('Cell Type')
set(gca, 'xtick', [1,2,3,4])
set(gca, 'xticklabel', {'ON parasol', 'ON midget', 'OFF parasol', 'OFF midget'})
ylabel('Biphasic Index')
title('Biphasic Index')
ylim = get(gca,'ylim');
axis([0 5 ylim])

peak_ratio_ON = peak_ratio(:,2)./peak_ratio(:,1)    
figure; hist(peak_ratio_ON)
xlabel('ON Parasol Peak Ratio/ON Midget Peak Ratio')
ylabel('Number of datasets')
title('Peak Ratio ON Cells')

peak_ratio_OFF = peak_ratio(:,4)./peak_ratio(:,3)    
figure; hist(peak_ratio_OFF)
xlabel('OFF Parasol Peak Ratio/OFF Midget Peak Ratio')
ylabel('Number of datasets')
title('Peak Ratio OFF Cells')


peak_to_zero_cross_ON = peak_to_zero(:,2)./peak_to_zero(:,1)    
figure; hist(peak_to_zero_cross_ON)
xlabel('ON Parasol Peak to Zero Crossing/ON MidgetPeak to Zero Crossing')
ylabel('Number of datasets')
title('Peak to Zero Crossing ON Cells')

peak_to_zero_cross_OFF = peak_to_zero(:,4)./peak_to_zero(:,3)    
figure; hist(peak_to_zero_cross_OFF)
xlabel('OFF Parasol Peak to Zero Crossing/OFF Midget Peak to Zero Crossing')
ylabel('Number of datasets')
title('Peak to Zero Crossing OFF Cells')

figure
plot(zero_crossing(:,1), peak_to_zero(:,1), 'xk')
hold on
plot(zero_crossing(:,2), peak_to_zero(:,2), 'or')
legend('ON Parasol', 'ON Midget')
xlabel('Zero Crossing')
ylabel('Peak to Zero')
title('Compare ON parasols and midgets on 2 dimensions')


figure
plot(zero_crossing(:,1), biphasic_index(:,1), 'xk')
hold on
plot(zero_crossing(:,2), biphasic_index(:,2), 'or')
plot(zero_crossing(:,3), biphasic_index(:,3), 'sg')

legend('ON Parasol', 'ON Midget', 'OFF parasol')
xlabel('Zero Crossing')
ylabel('Biphasic Index')
title('Compare parasols and midgets on 2 dimensions')



figure
plot(zero_crossing(:,3), peak_to_zero(:,3), 'xk')
hold on
plot(zero_crossing(:,4), peak_to_zero(:,4), 'or')
legend('OFF Parasol', 'OFF Midget')
xlabel('Zero Crossing')
ylabel('Peak to Zero')
title('Compare OFF parasols and midgets on 2 dimensions')


figure
plot(trough_to_zero(:,1), peak_to_zero(:,1), 'xk')
hold on
plot(trough_to_zero(:,2), peak_to_zero(:,2), 'or')
legend('ON Parasol', 'ON Midget')
xlabel('Trough to Zero')
ylabel('Peak to Zero')
title('Compare ON parasols and midgets on 2 dimensions')

figure
plot(trough_to_zero(:,3), peak_to_zero(:,3), 'xk')
hold on
plot(trough_to_zero(:,4), peak_to_zero(:,4), 'or')
legend('OFF Parasol', 'OFF Midget')
xlabel('Trough to Zero')
ylabel('Peak to Zero')
title('Compare OFF parasols and midgets on 2 dimensions')



figure
plot(zero_crossing(:,3), peak_ratio(:,3), 'xk')
hold on
plot(zero_crossing(:,4), peak_ratio(:,4), 'or')
legend('OFF Parasol', 'OFF Midget')
xlabel('Zero Crossing')
ylabel('Peak Ratio')
title('Compare OFF parasols and midgets on 2 dimensions')


figure
plot(zero_crossing(:,1), peak_ratio(:,1), 'xk')
hold on
plot(zero_crossing(:,2), peak_ratio(:,2), 'or')
legend('ON Parasol', 'ON Midget')
xlabel('Zero Crossing')
ylabel('Peak Ratio')
title('Compare ON parasols and midgets on 2 dimensions')