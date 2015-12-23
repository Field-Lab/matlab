% view ouput of compute_average_movie


% PARAMETERS
scale_up_factor = 1;


the_movie = average_movie;

avg_movie = the_movie.frames;
frame_count = the_movie.frame_count;
avg_frame = the_movie.avg_frame;
frame_offset = the_movie.frame_offset;
range = the_movie.range;




% select movie
switch 2
    case 1 % spike triggered movie


        % normalize movie

        % duplicate average movie
        norm_movie = avg_movie;

        % compute average frame
        %avg_frame = mean(norm_movie,3);

        % normalize each frame
        for ff=1:length(frame_count)

            % divide by the average
            %norm_movie(:,:,ff) = norm_movie(:,:,ff) ./ avg_frame;

            % set median to 0
            norm_movie(:,:,ff) = norm_movie(:,:,ff) - median(reshape(norm_movie(:,:,ff),[],1));
        end

        display_movie = norm_movie;

    case 2 % subsampled movie

        display_movie = the_movie.sub_movie;
end




% normalize amplitudes to fit in [0 64] window
%frames_min = min(reshape(display_movie,[],1));
%frames_max = max(reshape(display_movie,[],1));
vals = reshape(display_movie,[],1);
frames_min = median(vals) - max(vals - median(vals));
frames_max = median(vals) + max(vals - median(vals));
display_movie = 64 * (display_movie - frames_min)/(frames_max - frames_min);

display_movie = scale_up_factor * (display_movie - 32) + 32;




start_index=1; index_min=1; index_max=size(display_movie,3);

% initialize ROI
new_roi = false;
roi_im = [];
roi_ims = cell(0);
roi_pts = struct;
roi_intensity = [];

ri_colors = 'rgbmc';

figure(4);clf;colormap gray
slider = make_loop_slider_getpts(start_index,index_min,index_max);
first_time = true;

while 1

    ff = round(get(slider,'Value'));

    % frame image
    roi_axes = subplot('position',[0.01 0.1 0.58 .85]);cla(roi_axes)

    % note xlim, ylim
    if first_time
        % if first time, set to be relevant region
        curr_lim = [1 size(avg_frame,2); 1 size(avg_frame,1)];
        first_time = false;
    else
        curr_lim = [xlim; ylim];
    end


    % plot image
    image(display_movie(:,:,ff))
    axis image; hold on

    % restore bounds
    xlim(curr_lim(1,:))
    ylim(curr_lim(2,:))


    % add title
    title(sprintf('average of frames beginning between %0.1f and %0.1f msec (%d images)',...
        (range(1)+ff-1)*frame_offset*1000,(range(1)+ff)*frame_offset*1000,frame_count(ff)))


    % histogram of values
    hist_axes = subplot('position',[0.625 .5 0.35 .45]);

    % if no ROI, just plot histogram for this frame
    if isempty(roi_im)
        axes(hist_axes)
        hist(reshape(display_movie(:,:,ff),[],1),50)
        xlim([0 64])
        title('values in this frame')

    else
        % otherwise, do lots of things


        % compute signal in ROI
        if new_roi
            ri = zeros(1,size(display_movie,3));
            for gg = 1:size(display_movie,3)
                ri(gg) = sum(sum(display_movie(:,:,gg).*roi_im));
            end
            ri = ri / sum(sum(roi_im));

            % save
            roi_intensity(size(roi_intensity,1)+1,:) = ri;
            roi_pts(size(roi_intensity,1)).x = roi_x;
            roi_pts(size(roi_intensity,1)).y = roi_y;
            roi_ims{size(roi_intensity,1)} = roi_im;

            new_roi = false;
        end

        % make axes
        plot_axes = subplot('position',[0.625 .1 0.35 .35]);

        % set up axes
        cla(plot_axes)
        hold(plot_axes,'on')
        cla(hist_axes)
        hold(hist_axes,'on')

        % get frame
        this_frame = display_movie(:,:,ff);

        % plot vertical bar showing depth
        plot([ff ff],[min(min(roi_intensity)) max(max(roi_intensity))],'-','Color',.8*[1 1 1],'Parent',plot_axes)

        % plot stuff
        for rr = 1:size(roi_intensity,1)

            % plot intensity
            plot(roi_intensity(rr,:),'.-','color',ri_colors(rr),'parent',plot_axes)

            % draw ROI
            plot(roi_pts(rr).x,roi_pts(rr).y,'color',ri_colors(rr),'parent',roi_axes)

            % histogram
            [n,x] = hist(reshape(this_frame(roi_ims{rr}),[],1),0:2:64);
            plot(x,n/max(n),'color',ri_colors(rr),'parent',hist_axes)

        end

        % titles
        title(plot_axes,'average value in the ROI')

        % plot histogram of intensity in the ROI
        axes(hist_axes)
        xlim([0 64])
        title('values in the ROI')
    end


    % end with axes where ROI will be drawn
    axes(roi_axes)

    uiwait;
end


