% get timecourse properties

function [rf, t_zc, t_p, t_t, bi_ind, fr, amp] = get_timecourse_prop(datarun, cell_id, run_opts)

indices = get_cell_indices(datarun, cell_id);
sta = datarun.stas.stas{indices};
% try
%     [threshold] = sig_stixels_threshold(sta, params.false_stixels);
% catch
threshold = 3.5;
% end

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', threshold);

num_spikes = length(datarun.spikes{indices});


if sum(full(sig_stixels(:)))/length(full(sig_stixels(:))) < 0.003
    rf = nan;
    t_zc = nan;
    t_p = nan;
    t_t = nan;
    bi_ind = nan;
    amp = nan;
else
    
    tc = time_course_from_sta(sta, sig_stixels);
    
    num_colors = size(tc,2);
    if num_colors == 1
        color_index = 1;
    else
        color_index = 2;
    end
    
    % outlier_idx = abs(timecourse_green(:,27) - median(timecourse_green(:,27))) > 3*std(timecourse_green(:,27)); % Find outlier idx
    
    % timecourses_good = timecourse_green(~outlier_idx,:);
    % avg_timecourse = mean(timecourses_good,1)';
    % if abs(min(avg_timecourse)) > abs(max(avg_timecourse))
    %     avg_timecourse_orient  = -avg_timecourse;
    % else
    %     avg_timecourse_orient = avg_timecourse;
    % end
    
    time = -(run_opts.num_frames-run_opts.frames_past_zero-1)*run_opts.refresh:run_opts.refresh:run_opts.refresh + run_opts.frames_past_zero - 1;
    
    
    [max_ , max_ind] = max(tc(:,color_index));
    [min_ , min_ind] = min(tc(:,color_index));

    if max_ind > min_ind
        peak_ind = max_ind;
        trough_ind = min_ind;
        peak = max_;
        trough = min_;
        
    else
        peak_ind = min_ind;
        trough_ind = max_ind;
        peak = min_;
        trough = max_;
    end
    
    t_p = abs(peak_ind*run_opts.refresh + time(1)- run_opts.refresh);
    amp = tc(peak_ind,color_index);
    
    [t_ind,t0] = crossing(tc(:,color_index));
    
    less_than = find(t_ind < peak_ind == 1);
    if isempty(less_than)
        t_zc = nan;
    else
        
        ind = less_than(end);
        % [~,ind] = min(abs(t0 - run_opts.num_frames/2)); % find the frame that closest to the middle of num_frames
        t0 = t0(ind);
        t_zc = abs(t0*run_opts.refresh + time(1)- run_opts.refresh);
    end
    
%     if amp > 0
%     [trough , trough_ind] = min(tc(:,color_index));
%     else
%         [trough , trough_ind] = max(tc(:,color_index));
%     end
    
    t_t = abs(trough_ind*run_opts.refresh + time(1)- run_opts.refresh);
    
    bi_ind = abs(amp)/abs(trough);
    

    % gmean = nan(length(indicies),1);
    % for i = 1:length(indicies)
    gmean = geomean([datarun.stas.fits{indices}.sd(1)*datarun.stimulus.stixel_width*5.5,datarun.stas.fits{indices}.sd(2)*datarun.stimulus.stixel_width*5.5]); % 5.5 um/pixel
    %
    %
    % end
    
    
    rf = gmean*2; % want the diameter instead of the radius
end

fr = num_spikes/datarun.duration;