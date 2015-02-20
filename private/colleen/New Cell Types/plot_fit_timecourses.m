function plot_fit_timecourses(sta, fit_indices, fixed_indices, fit_params, fixed_params)

all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);

temp_stix = significant_stixels(sta, 'time', 'std', 'select', 'thresh', 'thresh', 2.0, 'robust_std_method', 3); %changed from 3.5 to 3.0
biggestBlob = ExtractNLargestBlobs(full(temp_stix), 1);
temp_stix = biggestBlob;
fit_tc = time_course_from_sta(sta_fit, temp_stix);
norm_factor = max(abs(reshape(fit_tc, 1, [])));
fit_tc = fit_tc ./ norm_factor;
if size(sta_fit, 3) == 3
    plot(fit_tc(:,1), 'r')
    hold on
    plot(fit_tc(:,2), 'g')
    plot(fit_tc(:,3), 'b')
elseif size(sta_fit, 3) == 1
    plot(fit_tc, '--k')
    hold on
else
    error('dimensions of sta color is not recognized')
end