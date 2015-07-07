function h = plot_all_timecourses(datarun, cell_specification, params, run_opts)

[cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);

% Get interval for plotting timecourse
datarun = load_java_movie(datarun, run_opts.movie_spec);
interval = datarun.stimulus.interval;

% get cell type name
cell_type_number = find_cell_types(datarun, cell_specification);
if cell_type_number(1) > 0
    cell_type = datarun.cell_types{cell_type_number(1)}.name;
else
    cell_type = 'unknown type';
end

% make title
figure
for i = 1:length(cell_indices)
    
cell_id=get_cell_ids(datarun,cell_specification);

% plot time_course, if it exists
time_course = datarun.vision.timecourses(cell_indices(i));
time_course = [time_course.r, time_course.g, time_course.b];
if ~isempty(time_course)
    if params.normalize
        time_course = time_course ./ norm(time_course);
    end
    h = plot_time_course_(time_course)
    hold on
%     xlabel('frame number')
%     ylabel('contrast')
end
end
x = 0:8.33*interval:8.33*interval*(size(datarun.vision.timecourses(1).r,1)-1);
set(gca, 'xtick',1:size(datarun.vision.timecourses(1).r,1))

set(gca, 'xticklabel',round(fliplr(x(1:end))))
title(sprintf('%s',cell_type_name))
% close(gcf)

