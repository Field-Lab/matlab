% Find midget and parasol zero crossing and ratio of peaks

% load data


file_name  = '2010-09-24-0/data001-nwpca/data001-nwpca';
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];

%% Load Data1
opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',0);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames =1:30;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

[ON_parasol] = get_cell_indices(datarun, 'ON midget');

timecourses = zeros(length(ON_parasol), 30);
for i = 1:length(ON_parasol)
    timecourses(i,:) = datarun.vision.timecourses(ON_parasol(i)).g;
end


outlier_idx = abs(timecourses(:,27) - median(timecourses(:,27))) > 3*std(timecourses(:,27)); % Find outlier idx
timecourses_good = timecourses(~outlier_idx,:);
avg_timecourse = mean(timecourses_good,1)';



time = [-29*datarun.stimulus.interval*8.33:datarun.stimulus.interval*8.33:0]';

figure
plot(time,avg_timecourse)
hold on
plot(time, zeros(size(time)), 'k')

[ind,t0] = crossing(avg_timecourse);
zero_crossing = (30-t0)*datarun.stimulus.interval*8.33

peak_ratio = abs(max(avg_timecourse))/abs(min(avg_timecourse))