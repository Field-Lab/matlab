function plot_sta_(sta, varargin)
% plot_sta_     plot an STA in an interactive slider window
%
% usage:  plot_sta_(sta, varargin)
%
% arguments:      sta - YxXxCxT matrix
%            varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% figure          	0            	figure or axes where to plot it. 0 for new figure, -1 in current axes
% prefix            ''              text for the title
% overlay           []              Nx2 matrix of points to overlay
%                                       usually the array outline
%
%
% 2009-06  gauthier
% 2010-01  phli
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);%, @(x)~isempty(x)&&round(x)==x);
p.addParamValue('prefix', '', @ischar);
p.addParamValue('overlay', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.figure);

% choose first frame to show
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));

% normalize STA color
sta = norm_image(sta);

% create slider control
ha = make_loop_slider_list(start_index, 1, size(sta, 4), {@slider_plot, sta, params});

% plot once before any clicks
slider_plot(ha, [], sta, params);


function slider_plot(handle, event, sta, params) %#ok<INUSL>
% display one frame of the STA

% get the slider position
ss = round(get(handle,'Value'));
cla;

% plot spatial sensitivity
image(sta(:,:,:,ss));
axis image
%set(gca,'xtick',[],'ytick',[])

if ~isempty(params.overlay)
    hold on
    plot(params.overlay(:,1),params.overlay(:,2),'k')
end

% title
title(sprintf('%sframe %d of %d (%d)',params.prefix,ss,size(sta,4),ss-size(sta,4)))
