function datarun = compute_monitor_to_array_transformation(datarun, varargin)
% compute_monitor_to_array_transformation     register array coordinates and monitor coordinates
%
% NOTE: for graphical documentation, see "matlab-standard/documentation/compute monitor to array transformation
%
%
%   MODE 1: if datarun.piece.array exists and is nonempty, these numbers are used to generate
%           an approximate transformation
%
%   MODE 2: otherwise, the user must supply photographic mapping images in the struct
%           datarun.piece.photographic_mapping.images (see below).  then the user will be prompted to select
%           alignment points for mapping the photographic mapping stimulus onto "pm#", optionally
%           "array" onto "array_edges", and finally "array_edges" onto a picture of the array.
%           The optional step is necessary because sometimes the image of the electrodes is zoomed in
%           and it's not clear what part of the array is being shown.  Note that "pm#" and the first array image
%           are assumed to be in the same coordinate frame, i.e. images of the same region focused on the stimulus
%           and electrodes, respectively.  If the user has already selected alignment points, these
%           will be stored in datarun and used as starting conditions when this function is called again.
%           
%           NOTE: AT LEAST FOUR ALIGNMENT POINTS MUST BE SELECTED OR THE CODE WILL CRASH!!!
%
%
%
% usage:  datarun = compute_monitor_to_array_transformation(datarun, <params>)
%
% arguments:     datarun - datarun struct with the field:
%
%                           datarun.piece.array - array "corners" from manual mapping
%
%                        OR with fields containing photographic mapping images:
%
%                           datarun.piece.photographic_mapping.images.pm#  [# = 2, 10, or 32]
%                               an image showing the photographic mapping stimulus
%                           datarun.piece.photographic_mapping.images.array
%                               optional, an image showing the electrodes in the same coordinates as "pm#"
%                           datarun.piece.photographic_mapping.images.array_edges
%                               an image in which electrodes numbers can be clearly identified,
%                               usually because a corner the array is visible.  if "array" is not
%                               provided, this is assumed to be in the same coordinates as "pm#"
%
%                       OR with fields specifying the paths to photographic mapping images:
%
%                           datarun.piece.photographic_mapping_pm#  [# = 2, 10, or 32]
%                           datarun.piece.photographic_mapping_array (optional)
%                           datarun.piece.photographic_mapping_array_edges
%
%                           the first time this function is called, the images will be loaded and stored in datarun
%
%
%               <params> - struct or list of optional parameters (see below)
%
%
% outputs:     datarun - datarun struct with fields:
%
%           datarun.piece.corners                                - Nx2 matrix, "corners" of the array in monitor coordinates
%                        .T_monitor_to_array                     - matlab tform struct
%                        .photographic_mapping.T_base_to_array   - matlab tform struct (mode 2 only)
%                                             .T_base_to_monitor - matlab tform struct (mode 2 only)
%                                             .points_*          - points used for mapping between various images (mode 2 only)
%
%           datarun.ei.position                                  - Mx2 matrix, x-y electrode locations in array coordinates
%                     .position_monitor                          - Mx2 matrix, x-y electrode locations in monitor coordinates
%                     .position_sta                              - Mx2 matrix, x-y electrode locations in STA coordinates (if white noise run)
%
%
%
% optional params, their default values, and what they specify:
%
% fig               []          	figure to plot in.  if 0, make new figure. if empty, don't plot.
% source            []              what to use to identify the transformation
%                                       'clicks' - use clicked array edges from standard mapping (mode 1)
%                                       'images' - use photographic mapping images with alignment points selected by the user (mode 2)
% which_pm        	'pm32'          which image to use for alignment (mode 2 only)
% reselect          true            if alignment points already exist, should the user be prompted to alter/approve them?
%                                       if true, a GUI will appear for the user to select alignment points
%
%
%
% examples:
%
%   % align pm32 to get coarse alignment
%   compute_monitor_to_array_transformation(datarun,'which_pm','pm32');
%
%   % align pm2 to get points at a finer scale
%   compute_monitor_to_array_transformation(datarun,'which_pm','pm2');
%
%
% 2010-01  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% note valid names for alignment images
valid_image_names = load_alignment_image;

% specify list of optional parameters
p.addParamValue('fig', []);
p.addParamValue('source',[],       @(x) any(strcmpi(x, {'images' 'clicks'})));
p.addParamValue('which_pm','pm32', @(x) any(strcmpi(x, valid_image_names)));
p.addParamValue('reselect', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% IDENTIFY AND VERIFY THE SOURCE OF ALIGNMENT INFORMATION


% identify what is present in datarun

% look for photographic mapping images
has_images = false;
if isfield(datarun,'piece') && isfield(datarun.piece,'photographic_mapping') && ...
        isfield(datarun.piece.photographic_mapping,'images') && ...
        any(isfield(datarun.piece.photographic_mapping.images,valid_image_names)) && ...
        isfield(datarun.piece.photographic_mapping.images,'array_edges')
    has_images = true;
elseif isfield(datarun, 'stacks') && size(datarun.stacks, 1) >= 4 && any(datarun.stacks(4,3:end))
    has_images = true;
end

% look for clicked mapping information
if isfield(datarun,'piece') && isfield(datarun.piece,'array') && ~isempty(datarun.piece.array)
    has_clicks = true;
else has_clicks = false;
end


% if a source was specified, verify the required data exist
if ~isempty(params.source)
    switch params.source
        case 'images'
            if ~has_images;
                fprintf('\n\n\n')
                disp(valid_image_names)
                error('datarun.piece.photographic_mapping must contain a field named "array" and a field with one of the names above.')
            end
        case 'clicks'
            if ~has_clicks;error('datarun.piece.array must exist and be non-empty.')
            end

    end

    % if no source was specified, check for what's available
else
    if has_images   % use photographs, if available
        params.source = 'images';

    elseif has_clicks    % otherwise use map information, if available
        params.source = 'clicks';

    else   % if neither are available, give error
        error('datarun struct does not contain the data required to perform image alignment.  see help for requirements.')
    end
end




% LOAD INFORMATION NEEDED IN ALL MODES


% get transformation from monitor to camera coordinates using rig/optical path information
T_monitor_to_camera = identify_monitor_to_camera_transform(datarun);

% get transformation from array image to camera coordinates using rig/optical path information
T_array_image_to_camera = identify_array_image_to_camera_transform(datarun);

% get array info
array_info = load_array_info(datarun,2);






% do the alignment already
if strcmp(params.source,'clicks')

    %%%%%%%%%%%%%%%%%%%%%%%%%% MODE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use clicked array edges to generate transformations


    % get where array corners are expected to be
    ideal_corners = array_info.corners;
    % transform them to camera coordinates
    ideal_corners_cam = tformfwd(maketform('composite',T_array_image_to_camera,array_info.T_array_to_array_image),ideal_corners);
    % center them
    ideal_corners_cam = ideal_corners_cam - repmat(mean(ideal_corners_cam),size(ideal_corners_cam,1),1);
    % set norm to 1
    ideal_corners_cam = ideal_corners_cam / norm(ideal_corners_cam);


    % get clicked array corners
    map_points = datarun.piece.array;
    % account for offset (obvius starts at [0,0], matlab at [1,1])
    map_points = map_points + 1;
    % transform them to camera coordinates
    map_points_cam = tformfwd(T_monitor_to_camera,map_points);
    % center them
    map_points_cam = map_points_cam - repmat(mean(map_points_cam),size(map_points_cam,1),1);
    % set norm to 1
    map_points_cam = map_points_cam / norm(map_points_cam);


    % visually compare the clicked points to the ideal corners
    %figure(1);clf;plot(ideal_corners_cam(:,1),ideal_corners_cam(:,2),'.k',map_points_cam(:,1),map_points_cam(:,2),'.r')


    % verify there are the correct number of mapped points
    if ~all(size(map_points) == [size(ideal_corners,1) 2])
        error('datarun.piece.array has the wrong size.  Should be [%d %d], is [%d %d]',...
            size(ideal_corners,1),2,size(map_points,1),size(map_points,2))
    end


    % identify correspondence between clicked corners and array corners
    % this calls the ugly but functional "TPS-RPM" package

    % get the [undocumented] matrix m
    [junk,junk,m] = cMIX(map_points_cam,ideal_corners_cam, 1, 0.1, 100,false,'icp3');

    % interpret it
    [junk,jj]=max(m,[],2);
    

    % compute transformation and save in datarun
    T_monitor_to_array = cp2tform(map_points,ideal_corners(jj,:),'projective');
    datarun.piece.T_monitor_to_array = T_monitor_to_array;
    
    

    % plot
    if 0
        figure(1);clf
        % get original mapping stimulus
        pmso = load_alignment_image('pm','pm32');
        % plot it
        imagesc(pmso.im,'xdata',pmso.xdata,'ydata',pmso.ydata);axis image;hold on;colormap gray
        % transform electrodes to monitor coordinates
        epm=tforminv(T_monitor_to_array,array_info.positions);
        % plot them
        plot(epm(:,1),epm(:,2),'.')
        % with numbered labels
        for ee=1:size(epm,1);
            text(epm(ee,1),epm(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',12,'HorizontalAlignment','Center','VerticalAlignment','Bottom');
        end
    end




elseif strcmp(params.source, 'images')
    datarun = load_pm_images(datarun);
    pm_images = get_pm_images(datarun);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% MODE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prompt user to select alignment points



    % LOAD PHOTOGRAPHIC MAPPING STIMULUS


    % load photographic mapping image in monitor coordinates
    pm_stim_mon = load_alignment_image('pm',params.which_pm);

    % apply transformation to image and x/ydata
    [pm_stim.im,pm_stim.xdata,pm_stim.ydata] = imtransform(pm_stim_mon.im,T_monitor_to_camera,'udata',pm_stim_mon.xdata,'vdata',pm_stim_mon.ydata);

    % compute transformation from monitor to pm_stim.im
    T_monitor_to_pm_stim = maketform('composite',maketform('affine',...
        [pm_stim.xdata(1) pm_stim.ydata(1); pm_stim.xdata(1) pm_stim.ydata(2); pm_stim.xdata(2) pm_stim.ydata(2)],...
        [1 1; 1 size(pm_stim.im,1); size(pm_stim.im,2) size(pm_stim.im,1)]),...
        T_monitor_to_camera);



    % LOAD IMAGE OF THE ARRAY

    % transform array image to camera coordinates
    array.im = imtransform(array_info.image,T_array_image_to_camera);
    
    % normalize images to be double with a peak of 1.0
    im_names=fieldnames(pm_images);
    for ii=1:length(im_names)
        pm_images.(im_names{ii}) = double(pm_images.(im_names{ii})); % added line
        pm_images.(im_names{ii})=pm_images.(im_names{ii})/max(pm_images.(im_names{ii})(:));
    end

    
    
    % REGISTER THE PHOTOGRAPHIC MAPPING STIMULUS WITH THE FIRST PM IMAGE


    % prompt user to select points

    % determine whether points already exist, and if so start with them
    if isfield(datarun.piece,'photographic_mapping') && isfield(datarun.piece.photographic_mapping,'points_camera_to_base') && ...
            isfield(datarun.piece.photographic_mapping.points_camera_to_base,'in') && ...
            isfield(datarun.piece.photographic_mapping.points_camera_to_base,'out')

        % only show the GUI if desired
        if params.reselect
            prompt_to_select_points_cmtaf
            
            [xy_pm_stim, xy_pm] = cpselect(pm_stim.im,pm_images.(sprintf('photographic_mapping_%s',params.which_pm)),...
                datarun.piece.photographic_mapping.points_camera_to_base.in,...
                datarun.piece.photographic_mapping.points_camera_to_base.out,...
                'wait',true);
        else
            xy_pm_stim = datarun.piece.photographic_mapping.points_camera_to_base.in;
            xy_pm = datarun.piece.photographic_mapping.points_camera_to_base.out;
        end
        
    else
        prompt_to_select_points_cmtaf
        
        [xy_pm_stim, xy_pm] = cpselect(pm_stim.im,pm_images.(sprintf('photographic_mapping_%s',params.which_pm)),'wait',true);
    end

    % store points
    datarun.piece.photographic_mapping.points_camera_to_base.in = xy_pm_stim;
    datarun.piece.photographic_mapping.points_camera_to_base.out = xy_pm;

    % compute and store transformation
    datarun.piece.photographic_mapping.T_camera_to_base = cp2tform(xy_pm_stim,xy_pm,'projective');




    % REGISTER "ARRAY" IMAGE WITH "ARRAY EDGES" IMAGE, IF NEEDED

    if isfield(pm_images,'array')

        % prompt user to select points

        % determine whether points already exist, and if so start with them
        if isfield(datarun.piece.photographic_mapping,'points_base_to_array') && ...
                isfield(datarun.piece.photographic_mapping.points_base_to_array,'in') && ...
                isfield(datarun.piece.photographic_mapping.points_base_to_array,'out')

            % only show the GUI if desired
            if params.reselect
                prompt_to_select_points_cmtaf
                
                [xy_array, xy_array_edges] = cpselect(pm_images.array,pm_images.array_edges,...
                    datarun.piece.photographic_mapping.points_base_to_array(1).in,...
                    datarun.piece.photographic_mapping.points_base_to_array(1).out,...
                    'wait',true);
            else
                xy_array = datarun.piece.photographic_mapping.points_base_to_array(1).in;
                xy_array_edges = datarun.piece.photographic_mapping.points_base_to_array(1).out;
            end
            
        else
            prompt_to_select_points_cmtaf
            
            [xy_array, xy_array_edges] = cpselect(pm_images.array,pm_images.array_edges,'wait',true);
        end

        % store points
        datarun.piece.photographic_mapping.points_base_to_array(1).in = xy_array;
        datarun.piece.photographic_mapping.points_base_to_array(1).out = xy_array_edges;

        % compute and store transformation
        datarun.piece.photographic_mapping.T_base_to_array(1) = cp2tform(xy_array,xy_array_edges,'projective');

    else
        % if pm.images.array doesn't exist, then just make this the identity transform
        datarun.piece.photographic_mapping.T_base_to_array(1) = maketform('affine',[0 1;0 0;1 0],[0 1;0 0;1 0]);

    end



    % REGISTER "ARRAY EDGES" IMAGE WITH THE ARRAY

    % prompt user to select points

    % determine whether points already exist, and if so start with them
    if isfield(datarun.piece.photographic_mapping,'points_base_to_array') && ...
            isfield(datarun.piece.photographic_mapping.points_base_to_array,'in') && ...
            isfield(datarun.piece.photographic_mapping.points_base_to_array,'out')

        % only show the GUI if desired
        if params.reselect
            prompt_to_select_points_cmtaf
            
            [xy_array_edges, xy_array] = cpselect(pm_images.photographic_mapping_array_edges,array.im,...
                datarun.piece.photographic_mapping.points_base_to_array(2).in,...
                datarun.piece.photographic_mapping.points_base_to_array(2).out,...
                'wait',true);
        else
            xy_array_edges = datarun.piece.photographic_mapping.points_base_to_array(2).in;
            xy_array = datarun.piece.photographic_mapping.points_base_to_array(2).out;
        end

    else
        prompt_to_select_points_cmtaf
            
        %[xy_array_edges, xy_array] = cpselect(pm_images.array_edges,array.im,'wait',true);
        [xy_array_edges, xy_array] = cpselect(pm_images.photographic_mapping_array_edges,array.im,'wait',true);
    end

    % store points
    datarun.piece.photographic_mapping.points_base_to_array(2).in = xy_array_edges;
    datarun.piece.photographic_mapping.points_base_to_array(2).out = xy_array;

    % compute and store transformation
    datarun.piece.photographic_mapping.T_base_to_array(2) = cp2tform(xy_array_edges,xy_array,'projective');




    % GENERATE TRANSFORMATIONS


    % from monitor coordinates to array coordinates
    T_monitor_to_array = maketform('composite',...
        fliptform(array_info.T_array_to_array_image),...
        fliptform(T_array_image_to_camera),...
        datarun.piece.photographic_mapping.T_base_to_array(2),...
        datarun.piece.photographic_mapping.T_base_to_array(1),...
        datarun.piece.photographic_mapping.T_camera_to_base,...
        T_monitor_to_pm_stim);

    % from base to array
    T_base_to_array = maketform('composite',...
        fliptform(array_info.T_array_to_array_image),...
        fliptform(T_array_image_to_camera),...
        datarun.piece.photographic_mapping.T_base_to_array(2),...
        datarun.piece.photographic_mapping.T_base_to_array(1));

    % from base to monitor
    T_base_to_monitor = maketform('composite',...
        fliptform(T_monitor_to_pm_stim),...
        fliptform(datarun.piece.photographic_mapping.T_camera_to_base));




    % PUT TRANSFORMATIONS INTO DATARUN

    datarun.piece.T_monitor_to_array = T_monitor_to_array;
    datarun.piece.photographic_mapping.T_base_to_array = T_base_to_array;
    datarun.piece.photographic_mapping.T_base_to_monitor = T_base_to_monitor;




    % PLOT, IF DESIRED

    % set up axes
    plot_axes = set_up_fig_or_axes(params.fig);

    % decide to plot or not
    if ~isempty(plot_axes)
        % get original mapping stimulus
        pmso = load_alignment_image('pm',params.which_pm);
        % plot it
        imagesc(pmso.im,'xdata',pmso.xdata,'ydata',pmso.ydata);axis image;hold on;colormap gray
        % transform electrodes to monitor coordinates
        epm=tforminv(T_monitor_to_array,array_info.positions);
        % plot them
        plot(epm(:,1),epm(:,2),'.')
        % with numbered labels
        for ee=1:size(epm,1);
            text(epm(ee,1),epm(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',12,'HorizontalAlignment','Center','VerticalAlignment','Bottom');
        end
    end


end



% COMPUTE AND SAVE ARRAY CORNERS USING THE TRANSFORMATION


% get array corners
corners = tforminv(T_monitor_to_array,array_info.corners);
% save in datarun
datarun.piece.corners = corners([1:end 1],:);
datarun.piece.corner_electrodes = array_info.corner_electrodes;



% COMPUTE AND SAVE ELECTRODE POSITIONS IN MONITOR AND STA COORDINATES

% array coordinates
datarun.ei.position = array_info.positions;
% monitor coordinates
datarun.ei.position_monitor = tforminv(T_monitor_to_array,array_info.positions);
% STA coordinates
if isfield(datarun,'stimulus') && isfield(datarun.stimulus,'stixel_width')
    T_monitor_to_STA = fliptform(coordinate_transform(datarun, 'monitor'));
    datarun.ei.position_sta = tformfwd(maketform('composite',[T_monitor_to_STA fliptform(T_monitor_to_array)]),array_info.positions);
end



function prompt_to_select_points_cmtaf

fprintf('Click to select alignment points.  YOU MUST SELECT AT LEAST FOUR POINTS!\n')


function datarun = load_pm_images(datarun)
datarun = copy_images_to_stacks(datarun);
for i = 1:length(datarun.stacks(4,:))
    datarun.stacks{4,i} = load_slices(datarun.stacks{4,i});
end


function pm_images = get_pm_images(datarun)
pm_images = struct();
for i = 1:length(datarun.stacks(4,:))
    if isempty(datarun.stacks{4,i})
        continue
    end
    stack = datarun.stacks{4,i};
    pm_images.(stack.name) = get_slice(stack, 1);
end