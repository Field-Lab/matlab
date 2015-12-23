% identify stratification of a given cell
% uses variable "stacks" to find best stack for this cell


% SELECT CELL

axon_id = 80;
roi_radius = 500;
plot_fig = 12;




% get soma center point
soma_center = axons{axon_id}(1,:);


% FIND THE BEST STACK

best_stack = identify_relevant_stacks(stacks,soma_center);
if isempty(best_stack)
    error('axon id %d is not contained in any of the stacks.  execute verify_fixed_image_alignment to check alignment.',axon_id)
else
    stack = stacks(best_stack(1));
end

stack = stacks(4);

% LOAD IMAGES
if 1
    
    % get location in stack coordinates
    roi_center = round(tforminv(stack.tform,soma_center(1,1),soma_center(1,2)));
    
    % set up ROI
    
    % get image info
    info = imfinfo(stack.im{1});
    % select part in the ROI
    roi_x = intersect(roi_center(1) + (-roi_radius:roi_radius),1:info.Width);
    roi_y = intersect(roi_center(2) + (-roi_radius:roi_radius),1:info.Height);
    % ensure it's ok
    if isempty(roi_x) || isempty(roi_y)
        fprintf('\n')
        error('no intersection with image\n    roi center = (%0.1f,%0.1f), image height = %d, width = %d',...
            roi_center(1),roi_center(2),size(im,1),size(im,2))
    end

    % load each image, and select roi
    clear ims
    for ii = 1:length(stack.im)
        fprintf('.')
        % load image in the ROI
        ims(:,:,:,ii) = imread(stack.im{ii},'PixelRegion', {[roi_y(1) roi_y(end)], [roi_x(1) roi_x(end)]});
    end
    fprintf('\n')
end



% PLOT STACK
if 1
    % set up figure
    figure(plot_fig);clf;
    start_index=round(size(ims,4)/2); index_min=1; index_max=size(ims,4);
    
    % initialize points
    pts = [];
    
    % transform axon to these coordinates
    axon_here = tforminv(stack.tform,axons{axon_id});
    
    % make slider
    slider = make_loop_slider_getpts(start_index,index_min,index_max);
    
    
    keep_going = 1;
    while keep_going
        k = round(get(slider,'Value'));
        
        % plot image at this plane
        pts_axes = subplot('Position',[.05 .1 .45 .85]);
        cla(pts_axes)
        imagesc(ims(:,:,:,k),'xdata',[min(roi_x) max(roi_x)],'ydata',[min(roi_y) max(roi_y)]);axis image; hold on
        % add the axons
        plot(axon_here(:,1),axon_here(:,2),'r')
        
        % plot points, if they exist
        if ~isempty(pts)
            % plot points
            plot(pts(:,1),pts(:,2),'ro','MarkerSize',20,'Parent',pts_axes)
            % make transformation to convert points to coordinates of the images
            im_x = size(ims,2);im_y = size(ims,1);
            T = maketform('affine',[min(roi_x) min(roi_y);min(roi_x) max(roi_y);max(roi_x) min(roi_y)],[1 1;1 im_y;im_x 1]);
            
            % plot depth profiles of the stack at each point
            
            % set up axes
            plot_axes = subplot_axes(plot_fig,[.55 .1 .45 .85],0,.05,1,size(pts,1));
            % for each point...
            for pp=1:size(pts,1)
                % plot its label
                text(pts(pp,1),pts(pp,2),num2str(pp),'Parent',pts_axes,'Color','r','FontSize',15)
                % get depth profile
                dp = depth_profile(ims,distance_from_point(size(ims(:,:,1,1)),tformfwd(T,pts(pp,:)),'radius',11)<10,'norm','max,min');
                % plot vertical bar showing depth
                hold(plot_axes{pp},'on')
                cla(plot_axes{pp})
                plot([k k],[0 max(max(dp))],'-','Color',.8*[1 1 1],'Parent',plot_axes{pp})
                % plot depth profile
                plot(fliplr(dp),'-','Parent',plot_axes{pp})
            end
        end

        % if user selected points, add them
        if exist('x','var');
            pts = [pts; x y];
            clear x y
        end
        
        % add title
        title(sprintf('axon %d, stacks(%d)',axon_id,best_stack))
        
        uiwait;
    end
    
end


