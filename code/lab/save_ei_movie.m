function save_ei_movie(datarun,cell_id,save_path,varargin)
% save_ei_movie     Save a movie of an EI
%
% usage:  save_ei_movie(datarun,cell_id,save_path, <params>)
%
% arguments:     datarun - datarun struct
%                cell_id - cell id
%              save_path - string, path to file
%               varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% frames            []         	which frames to plot, e.g. 10:40 .  if empty, plot every frame (except those that are all zero)
% figure            0          	figure to plot in. if 0, make new figure. if -1, use current figure.
% size              [1000 1000]	2-length vector, size of the movie.  if empty, use figure size
% fps               15          frames per second
% <other parameters>            any other parameters are passed to plot_ei.
%                                   by default, plot_ei is called with the following:
%
%                                   cutoff      -1
%                                   pos_color	[1 0 0]
%                                   neg_color	[0 0 1]
%
%
% 2010-01  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% some provided arguments might be for plot_ei and should be ignored.
p.KeepUnmatched = true;

% specify list of optional parameters
p.addParamValue('frames', []);
p.addParamValue('figure', 0);
p.addParamValue('size', [1000 1000]);
p.addParamValue('fps', 15);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
set_up_fig_or_axes(params.figure);
fig_num = gcf;

% set figure size
if isempty(params.size)
    position = get(fig_num,'position');
else
    % get current position of figure
    orig_position = get(fig_num,'position');
    % update x position
    position(1) = orig_position(1);
    position(3) = params.size(1);
    % update y position
    position(4) = params.size(2);
    position(2) = orig_position(2) + ( orig_position(4) - params.size(2));
end

% if user did not specify frames, use all frames, except those that are all zeros
if isempty(params.frames)
    % get the ei
    ei = get_ei(datarun,cell_id);
    % identify nonzero frames
    params.frames = find(~all(ei==0,1));
end


% initialize movie
MakeQTMovie('start',save_path);
MakeQTMovie('framerate',params.fps);
MakeQTMovie('quality',1);



% go through each frame
for ff=params.frames
    % prepare figure
    figure(fig_num); 
    clf;set(gcf, 'color', 'white')
    % make axes
    axes('position',[0 0 1 1]);
    % plot ei
    plot_ei(datarun,cell_id,'frame_number',ff,'cutoff',-1,'neg_color',[0 0 1],'pos_color',[1 0 0],p.Unmatched);
    set(gca, 'Visible', 'off');
    % remove clutter
    set(gca,'XTick',[],'YTick',[]);
    % ensure figure size is correct
    set(fig_num,'position',position);
    % ensure figure is drawn before adding to the movie
    drawnow
    % add to the growing list of frames
    MakeQTMovie('addaxes');
end

% finalize file
MakeQTMovie('finish')



