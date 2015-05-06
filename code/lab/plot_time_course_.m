function h = plot_time_course_(time_course, varargin)
% plot_time_course_     plot time course
%
% usage:  h = plot_time_course_(time_course, varargin)
%
% arguments:     time_course - TxC matrix, where C = 3 or 1
%                               T = # timesteps
%                               C = # colors
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     h - axis handle
%
%
% optional params, their default values, and what they specify:
%
% foa           	-1          figure or axes to plot in. if 0, make new figure. if -1, plot in current axes
% colors            'k'         for BW timecourse
%                   ['rgb']'    for RGB timecourse
%                               Cx3 double matrix or Cx1 string specifying plot color for each color
%
% exmaple:
%  
%     figure(1);clf;plot_time_course_(datarun.stas.time_courses{1})
%
%
%
% 2009-06 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', -1, @(x)~isempty(x));
switch size(time_course,2)
    case 1 % BW
        p.addParamValue('colors', 'k');
    case 3 % RGB
        p.addParamValue('colors', ['rgb']');
    otherwise % ??
        p.addParamValue('colors', rand(size(time_course,2),3));
end

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa, 0);

% plot each color
hold on
plot(plot_axes,[1 size(time_course,1)],[0 0],'color',.8*[1 1 1])
for cc=1:size(time_course,2)
    plot(plot_axes,(1:size(time_course,1))',time_course(:,cc),'-','Color',params.colors(cc,:))
end


if nargout > 0
    h = plot_axes;
end

