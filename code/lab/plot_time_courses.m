function plot_time_courses(datarun, cell_spec, varargin)
% plot_time_courses     plot time courses of several cells
%
% usage:  plot_time_courses(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%                  
%            varargin - struct or list of optional parameters (see below)
%
% cannot be used
%
% optional parameters, their default values, and what they specify:
%
%
% figure       	-1     	figure to plot in. if 0, make new figure. if -1, plot in current figure
% all           false   if false plots with slider, if 1 superimposes all plots
% bw            false   if true time courses are b&w, 
%                           
%
% 2009-06  gauthier
% 2012-08  sneha modified to include superimposing of time course plots
% 2013-03  sneha modified to plot normalized time courses when superimposed - normalized wrt norm
% 2014-08 GDF made changes to make all and bw inputs optional.


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure',0);
p.addParamValue('all', false, @islogical);
p.addParamValue('bw', false, @islogical);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
[junk,plot_fig] = set_up_fig_or_axes(params.figure);

% get cell ids
cell_ids = get_cell_ids(datarun,cell_spec);

%Slider plots
if ~p.Results.all
    % initialize figure with slider
    slider = make_loop_slider_list(1,1,length(cell_ids));
    % plot time course for each clel
    while 1
        % get the slider position
        cc = round(get(slider,'Value'));
        cla;
        % plot it
        plot_time_course(datarun,cell_ids(cc),'figure',-1,'clear',false)
        uiwait
    end
end

%Superimposed plots
if p.Results.all
    %make title
    title('Time Courses of Cells');
    
    if ~p.Results.bw % If time courses have rgb columns they have to each be plotted separately and then superimposed
        for a = 1:length(cell_ids)
            cell_index = get_cell_indices(datarun,cell_ids(a));
            time_course = datarun.stas.time_courses{cell_index};
            if ~isempty(time_course)
                time_course = time_course./norm(time_course); %normalize time course
                plot_time_course_(time_course)
                xlabel('frame number')
                ylabel('contrast')
            end
        end
        
    elseif p.Results.bw
        timcs = get_time_courses_matrix(datarun, cell_ids);
        for b = 1:length(cell_ids)
            timcs(:,b) = timcs(:,b)./norm(timcs(:,b));             %normalize time course
        end 
        if ~isempty(timcs)
                plot_time_course_(timcs);
                xlabel('frame number');
                ylabel('contrast');
                set(findobj('Type','line'),'Color','b');
        end
    end   
end
