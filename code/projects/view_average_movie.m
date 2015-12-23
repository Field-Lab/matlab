% view ouput of compute_average_movie


% PARAMETERS


% choose movie

which_movie = 'sta';  % spike triggered average movie
%which_movie = 'submovie';  % sub movie (to check for motion)


% spatial filtering

filter_radius = 0;  % no filtering
%filter_radius = 2;


% type of normalization

normalization = 'median';  % set median of each frame to 0
%normalization = 'z score';  % use z score (pool values over time at a single pixel)


% scale factor for display
scale_up_factor = 10;





% rename variables
the_movie = average_movie;

avg_movie = the_movie.frames;
frame_count = the_movie.frame_count;
avg_frame = the_movie.avg_frame;
frame_offset = the_movie.frame_offset;
range = the_movie.range;




% select movie
switch which_movie
    case 'sta' % spike triggered movie

        % duplicate average movie
        display_movie = avg_movie;

        % normalize
        switch normalization
            case 'median'  % set median of each frame to 0
                
                for ff=1:length(frame_count)
                    display_movie(:,:,ff) = display_movie(:,:,ff) - median(reshape(display_movie(:,:,ff),[],1));
                end

            case 'z score'  % normalize each pixel to Z score based on time course
                
                mn = mean(display_movie,3);
                sd = std(display_movie,[],3);
                for ff=1:size(display_movie,3)
                    display_movie(:,:,ff) = (display_movie(:,:,ff)-mn)./sd;
                end
                
        end


    case 'submovie' % subsampled movie
        display_movie = the_movie.sub_movie;
end


% normalize amplitudes to fit in [0 64] with median 32
vals = reshape(display_movie,[],1);
vals = vals(~isnan(vals));
frames_min = median(vals) - max(vals - median(vals));
frames_max = median(vals) + max(vals - median(vals));
display_movie = 64 * (display_movie - frames_min)/(frames_max - frames_min);




% scale up
display_movie = scale_up_factor * (display_movie - 32) + 32;


start_index=1; index_min=1; index_max=size(display_movie,3);

% initialize ROI
new_roi = false;
roi_im = [];
roi_ims = cell(0);
roi_pts = struct;
roi_intensity = [];

ri_colors = 'rgbmc';

figure(9);clf;colormap gray
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


    % grab display frame
    disp_frame = display_movie(:,:,ff);
    
    % filter if desired
    if filter_radius
        disp_frame = imfilter(disp_frame,fspecial('gaussian',filter_radius*3*[1 1],filter_radius));
    end
    
    % plot image
    image(disp_frame)
    axis image; hold on

    % restore bounds
    xlim(curr_lim(1,:))
    ylim(curr_lim(2,:))


    % add title
    if strcmp(which_movie,'sta')
        title(sprintf('average of frames beginning between %0.1f and %0.1f msec (%d images)',...
            (range(1)+ff-1)*frame_offset*1000,(range(1)+ff)*frame_offset*1000,frame_count(ff)))
    end
    
    % note scale up factor
    xlabel(sprintf('image intensity scaled linearly by %d',scale_up_factor))
        

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
            switch 1
                case 1 % mean
                    for gg = 1:size(display_movie,3)
                        ri(gg) = sum(sum(display_movie(:,:,gg).*roi_im));
                    end
                    ri = ri / sum(sum(roi_im));
                case 2 % variance
                    for gg = 1:size(display_movie,3)
                        ri(gg) = var(reshape(display_movie(:,:,gg).*roi_im,[],1));
                    end
            end

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


