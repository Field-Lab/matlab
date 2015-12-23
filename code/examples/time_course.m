

% load datarun
if ~exist('datarun','var')
    datarun = load_data('2007-09-18-1/data000-mg/data000/data000');
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'verbose',1,'load_sta',{1});
    datarun = get_sta_summaries(datarun,{1});
end


% plot one timecourse
figure(1);clf;plot_time_course(datarun,18)

% compute and plot average timecourse
avg_time_course = average_time_course(datarun, {1});
figure(2);clf;plot_time_course_(avg_time_course);

% plot several timecourses
figure(3);clf;plot_time_courses(datarun,{1})
