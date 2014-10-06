function [profiles,rois] = select_stack_roi_(stack, varargin)
% select_stack_roi     plot an image stack, select ROIs, and plot depth profile for each color channel
%
% usage:  result = my_function(arg1, <params>)
%
% argument:    stacks - YxXxCxZ matrix of the image data
% 
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% figure            0                   figure to plot in. if 0, make new figure.
%                                           if -1, plot in current.
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
% date author
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);
p.addParamValue('foo','bar', @(x)any(strcmpi(x,{'bar','bore'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% INITIALIZE STUFF


% set up plot axes
set_up_fig_or_axes(params.figure);
plot_fig = gcf;

% initialize points
%pts = [];

% transform axon to these coordinates
%axon_here = tforminv(stack.tform,axons{axon_id});




% MAKE FIGURE

% set up figure
figure(plot_fig);clf;
start_index=round(size(stack,4)/2); index_min=1; index_max=size(stack,4);


% create slider control
ha = make_fig(start_index, index_min, index_max, {@plot_one_frame, stack, params});

% plot once before any clicks
plot_one_frame(ha, [], stack, params);



end




function plot_one_frame(handle, event, stack, params) %#ok<INUSL>
% display one image in the stack

% get the slider position
k = round(get(handle,'Value'));
cla;


% plot image at this plane
pts_axes = subplot('Position',[.05 .1 .45 .85]);
cla(pts_axes)
imagesc(stack(:,:,:,k))%,'xdata',[min(roi_x) max(roi_x)],'ydata',[min(roi_y) max(roi_y)]);
axis image; hold on
% add the axons
%plot(axon_here(:,1),axon_here(:,2),'r')

% get user data
ud = get(gcf,'UserData');

% plot roi, if it exists
if isfield(ud,'new_roi')

    ud.points
    
    plot(ud.points(:,1),ud.points(:,2),'r')

    set(gcf,'UserData')

    % plot points
    %         plot(pts(:,1),pts(:,2),'ro','MarkerSize',20,'Parent',pts_axes)
    %         % make transformation to convert points to coordinates of the images
    %         im_x = size(stack,2);im_y = size(stack,1);
    %         T = maketform('affine',[min(roi_x) min(roi_y);min(roi_x) max(roi_y);max(roi_x) min(roi_y)],[1 1;1 im_y;im_x 1]);
    %
    %         % plot depth profiles of the stack at each point
    %
    %         % set up axes
    %         plot_axes = subplot_axes(plot_fig,[.55 .1 .45 .85],0,.05,1,size(pts,1));
    %         % for each point...
    %         for pp=1:size(pts,1)
    %             % plot its label
    %             text(pts(pp,1),pts(pp,2),num2str(pp),'Parent',pts_axes,'Color','r','FontSize',15)
    %             % get depth profile
    %             dp = depth_profile(stack,distance_from_point(size(stack(:,:,1,1)),tformfwd(T,pts(pp,:)),'radius',11)<10,'norm','max,min');
    %             % plot vertical bar showing depth
    %             hold(plot_axes{pp},'on')
    %             cla(plot_axes{pp})
    %             plot([k k],[0 max(max(dp))],'-','Color',.8*[1 1 1],'Parent',plot_axes{pp})
    %             % plot depth profile
    %             plot(fliplr(dp),'-','Parent',plot_axes{pp})
    %         end
end

% if user selected points, add them
%if exist('x','var');
%    pts = [pts; x y];
%    clear x y
%end

% add title
%title(sprintf('axon %d, stacks(%d)',axon_id,best_stack))

%uiwait;

end





function slider = make_fig(start_index,index_min,index_max,callback)
% GUI slider for loop over plot
%

% button to select ROI
b1= uicontrol('Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[0.7 0 .1 .08], ...
    'String','add ROI',...
    'Callback','[roi_im, roi_x, roi_y] = roipoly;ud=get(gcf,''UserData'');ud.new_roi=struct(''roi'',roi_im,''points'',[roi_x roi_y]);set(gcf,''UserData'',ud); uiresume;');

% button to finish
b3= uicontrol('Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[0.9 0 .1 .08], ...
    'String','done',...
    'Callback','uiresume; keep_going = 0;');

% slider
slider = uicontrol(gcf,...
    'Style','slider',...
    'Min' ,index_min,'Max',index_max, ...
    'Units','normalized', ...
    'Position',[0,0,0.7,.04], ...
    'Value', start_index,...
    'SliderStep',1/(index_max-index_min) * [1 1],...
    'CallBack', callback);

set(gcf,'Toolbar','figure')

end


