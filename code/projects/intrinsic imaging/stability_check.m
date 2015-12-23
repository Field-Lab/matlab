function average_movie = stability_check(spikes,image_spec,varargin)
% compute_average_movie     compute average movie accompanying datarun
%
% usage:  average_movie = compute_average_movie(spikes,image_spec,<params>)
%
% arguments:     spikes - vector of spike times (sec)
%            image_spec - struct with fields:
%
%           	         image_spec.path                path to folder with images
%                        image_spec.data(N).name        name of image stack (relative to the path)
%                                          .index       which image in the stack
%                                          .time        time of image acquisition (sec)
%
%             <params> - struct or list of optional parameters (see below)
%
%
% outputs:  average_movie - struct with fields:
%
%          average_movie.frames      - YxXxT matrix
%                       .frame_count - T-length vector, # images per frame
%                       .avg_frame   - YxX matrix, average of all images (unweighted by spikes)
%                       .sub_movie   - subset of the original movie
%                       .frame_offset - value of input parameter
%                       .range       - value of input parameter
%                   
%
%
%
% optional params, their default values, and what they specify:
%
% verbose           false         	show output
% subsample         []              period of frames put into the subsampled movie
%                                       if empty, don't make subsampled movie
% frame_offset      0.01            size of time bin in which to combine images (sec)
% range             [-50 50]        range of time bins around the spike
%
%
% 2010-03  gauthier
%

 

% load one image
temp = imread([im_path im_data(1).name],im_data(1).index);

% initialize movie storage
avg_movie = zeros(size(temp,1),size(temp,2),num_images);


% select subset of images
image_index = round(1:(length(im_data)/num_images):length(im_data));

% go through each image
for jj = 1:num_images

    ii = image_index(jj);

    % get the image
    im = double(imread([im_path im_data(ii).name],im_data(ii).index));

    % get image time
    % trigger times are when the trigger begins for each image
    im_time = triggers(ii);

    % save the image
    avg_movie(:,:,jj) = im;

    % show progress
    if mod(jj,10)==0
        fprintf('loaded %d of %d in %0.2f sec\n',jj,num_images,toc)
    end
end



average_movie = avg_movie;




function [im,im_time] = load_one_image(image_spec,index)
% get image and its time

im = double(imread([image_spec.path image_spec.data(index).name],image_spec.data(index).index));

im_time = image_spec.data(index).time;




% % COMPUTE AVERAGE
% if 1
% 
% 
%     % PARAMETERS
% 
%     num_images = 100;
% 
%     datarun = load_data('/Analysis/Gauthier/2010-02-26-0/data000/data000');
% 
%     datarun = load_neurons(datarun);
% 
%     triggers = datarun.triggers([5:4004   4006:8005   8007:12006   12008:16007]);
%     %triggers = datarun.triggers(5:4004);
%     %triggers = triggers(randperm(length(triggers)));
% 
%     % image information
%     im_path = '/Analysis/Gauthier/2010-02-26-0/sequences/';
%     im_data = struct;
% 
%     [im_data(1:4000).name] = deal('s1.tif');
%     for ii=1:4000; im_data(ii).index = ii; end
% 
%     [im_data(4001:8000).name] = deal('s2.tif');
%     for ii=4001:8000; im_data(ii).index = ii-4000; end
% 
%     [im_data(8001:12000).name] = deal('s3.tif');
%     for ii=8001:12000; im_data(ii).index = ii-8000; end
% 
%     [im_data(12001:16000).name] = deal('s4.tif');
%     for ii=12001:16000; im_data(ii).index = ii-12000; end
% 
% 
% 
% 
%     % compute average
% 
%     tic
% 
% 
% 
%     % load one image
%     temp = imread([im_path im_data(1).name],im_data(1).index);
% 
%     % initialize movie storage
%     avg_movie = zeros(size(temp,1),size(temp,2),num_images);
% 
% 
%     % select subset of images
%     image_index = round(1:(length(im_data)/num_images):length(im_data));
% 
%     % go through each image
%     for jj = 1:num_images
% 
%         ii = image_index(jj);
% 
%         % get the image
%         im = double(imread([im_path im_data(ii).name],im_data(ii).index));
% 
%         % get image time
%         % trigger times are when the trigger begins for each image
%         im_time = triggers(ii);
% 
%         % save the image
%         avg_movie(:,:,jj) = im;
% 
%         % show progress
%         if mod(jj,10)==0
%             fprintf('loaded %d of %d in %0.2f sec\n',jj,num_images,toc)
%         end
%     end
% 
% end
% 
% 
% 
% 
% % VIEW MOVIE
% if 1
% 
%     % normalize movie
% 
%     % duplicate average movie
%     norm_movie = avg_movie;
% 
%     % normalize amplitudes to fit in [0 64] window
%     frames_min = min(reshape(norm_movie,[],1));
%     frames_max = max(reshape(norm_movie,[],1));
%     norm_movie = 64 * (norm_movie - frames_min)/(frames_max - frames_min);
% 
% 
% 
%     start_index=1; index_min=1; index_max=size(norm_movie,3);
% 
%     % initialize ROI
%     new_roi = false;
%     roi_im = [];
%     roi_ims = cell(0);
%     roi_pts = struct;
%     roi_intensity = [];
% 
%     ri_colors = 'rgbmc';
% 
%     figure(3);clf;colormap gray
%     slider = make_loop_slider_getpts(start_index,index_min,index_max);
%     while 1
% 
%         ff = round(get(slider,'Value'));
% 
%         % frame image
%         roi_axes = subplot('position',[0.01 0.1 0.58 .85]);cla(roi_axes)
%         imagesc(norm_movie(:,:,ff));hold on
%         axis image
% 
% 
%         % histogram of values
%         hist_axes = subplot('position',[0.625 .5 0.35 .45]);
% 
%         % if no ROI, just plot histogram for this frame
%         if isempty(roi_im)
%             axes(hist_axes)
%             hist(reshape(norm_movie(:,:,ff),[],1),50)
%             xlim([0 64])
%             title('values in this frame')
% 
%         else
%             % otherwise, do lots of things
% 
% 
%             % compute signal in ROI
%             if new_roi
%                 ri = zeros(1,size(norm_movie,3));
%                 for gg = 1:size(norm_movie,3)
%                     ri(gg) = sum(sum(norm_movie(:,:,gg).*roi_im));
%                 end
%                 ri = ri / sum(sum(roi_im));
% 
%                 % save
%                 roi_intensity(size(roi_intensity,1)+1,:) = ri;
%                 roi_pts(size(roi_intensity,1)).x = roi_x;
%                 roi_pts(size(roi_intensity,1)).y = roi_y;
%                 roi_ims{size(roi_intensity,1)} = roi_im;
% 
%                 new_roi = false;
%             end
% 
%             % make axes
%             plot_axes = subplot('position',[0.625 .1 0.35 .35]);
% 
%             % set up axes
%             cla(plot_axes)
%             hold(plot_axes,'on')
%             cla(hist_axes)
%             hold(hist_axes,'on')
% 
%             % get frame
%             this_frame = norm_movie(:,:,ff);
% 
%             % plot vertical bar showing depth
%             plot([ff ff],[min(min(roi_intensity)) max(max(roi_intensity))],'-','Color',.8*[1 1 1],'Parent',plot_axes)
% 
%             % plot stuff
%             for rr = 1:size(roi_intensity,1)
% 
%                 % plot intensity
%                 plot(roi_intensity(rr,:),'.-','color',ri_colors(rr),'parent',plot_axes)
% 
%                 % draw ROI
%                 plot(roi_pts(rr).x,roi_pts(rr).y,'color',ri_colors(rr),'parent',roi_axes)
% 
%                 % histogram
%                 [n,x] = hist(reshape(this_frame(roi_ims{rr}),[],1),0:2:64);
%                 plot(x,n/max(n),'color',ri_colors(rr),'parent',hist_axes)
% 
%             end
% 
%             % titles
%             title(plot_axes,'average value in the ROI')
% 
%             % plot histogram of intensity in the ROI
%             axes(hist_axes)
%             xlim([0 64])
%             title('values in the ROI')
%         end
% 
% 
%         % end with axes where ROI will be drawn
%         axes(roi_axes)
% 
%         uiwait;
%     end
% 
% 
% end
% 
% 
% 
