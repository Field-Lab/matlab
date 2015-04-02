function fit_tc = plot_fit_timecourses(sta, fit_indices, fixed_indices, fit_params, fixed_params, sig_stixels, plot_raw)

all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);

temp_stix = significant_stixels(sta, 'time', 'std', 'select', 'thresh', 'thresh', 4.0); %changed from 3.5 to 3.0
% biggestBlob = ExtractNLargestBlobs(full(temp_stix), 1);
% temp_stix = biggestBlob;
% temp_stix = sig_stixels
% fit_tc = time_course_from_sta(sta_fit, temp_stix);
fit_tc = time_course_from_sta(sta_fit, temp_stix);

norm_factor = max(abs(reshape(fit_tc, 1, [])));
fit_tc = fit_tc ./ norm_factor;
if size(sta_fit, 3) == 3
    plot(linspace(1,size(sta,4),size(fit_tc,1)), fit_tc(:,1), 'r')
    hold on
    plot(linspace(1,size(sta,4),size(fit_tc,1)),fit_tc(:,2), 'g')
    plot(linspace(1,size(sta,4),size(fit_tc,1)),fit_tc(:,3), 'b')
elseif size(sta_fit, 3) == 1
    plot(linspace(1,size(sta,4),size(fit_tc,1)), fit_tc, '--k')
    hold on
else
    error('dimensions of sta color is not recognized')
end

if plot_raw
% Raw Data
tc = time_course_from_sta(sta, temp_stix);

norm_factor = max(abs(reshape(tc, 1, [])));
tc = tc ./ norm_factor;
if size(sta_fit, 3) == 3
    hold on 
    plot(linspace(1,size(sta,4),size(tc,1)), tc(:,1), ':r')
    hold on
    plot(linspace(1,size(sta,4),size(tc,1)),tc(:,2), ':g')
    plot(linspace(1,size(sta,4),size(tc,1)),tc(:,3), ':b')
elseif size(sta_fit, 3) == 1
    plot(linspace(1,size(sta,4),size(tc,1)), tc, '--k')
    hold on
else
    error('dimensions of sta color is not recognized')
end
end

