function plot_fig = plot_average_rf_profile(datarun,cell_types, varargin)
% plot_average_rf_profile     plot average rf profile for several cell types
%
% usage:  plot_fig = plot_average_rf_profile(datarun,cell_types, varargin)
%
% arguments:     datarun - datarun struct
%             cell_types - numbers of the cell types to plot, e.g. [1:4]
%               varargin - struct or list of optional parameters (see below)
%
% outputs:     plot_fig - figure where the results were plotted
%
%
% optional params, their default values, and what they specify:
%
% figure            []            	figure to plot in. if 0, make new figure.
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%       passed to 'get_rf_cone_profiles'
%
%         	center_type       	center point of the RF
%         	normalize       	how to normalize the amplitude of RF inputs
%           radius              radius in which to plot points
%
%
%       passed to 'curve_from_binning'
%
%           bin_edges           how to bin cone weights
%           average_y           how to combine cone weights in a bin
%
%
% gauthier 2008-12
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);

% params to pass to get_rf_cone_profiles
p.addParamValue('normalize', 'default value');
p.addParamValue('center_type', 'default value');
p.addParamValue('radius', 'default value');
%       to curve_from_binning
p.addParamValue('average_y', 'default value');
p.addParamValue('bin_edges', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
profiles_params = make_struct_to_pass(p.Results,{'normalize','normalize',...
    'center_type','center_type','radius','radius'});
curve_params = make_struct_to_pass(p.Results,{'bin_edges','bin_edges','average_y','average_y'});



% BODY OF THE FUNCTION


% make new figure, if desired
if params.figure==0;params.figure = figure;end

% set parameters for nice printing
set(gcf,'PaperUnits','centimeters')
xSize = 18; ySize=28;
xLeft = (22-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])

% generate subplot axes
plot_axes = subplot_axes(params.figure,[.07 .05 .93 .9],.13,.2,2,3);



if ~isempty(plot_axes)
    %axes(plot_axes)
    
    for cc=1:length(cell_types)
        
        cell_type = cell_types(cc);
    
        % get average RF profile
        [all_x, all_y] = get_rf_profiles(datarun,{cell_type},...
            'radius',max(params.bin_edges),profiles_params);
        
        % put each color in a different cell index
        for dd=1:size(all_y,2)
            all_y_{dd} = all_y(:,dd);
            all_x_{dd} = all_x(:);
        end
        
        % interpolate a curve
        [average_x, average_y] = curve_from_binning(all_x_,all_y_,curve_params);
        
        
        % plot the curve(s)
        
        % clear axes
        axes(plot_axes{cc});cla
        
        % add a zero line
        plot([min(params.bin_edges) max(params.bin_edges)],[0 0],'-','Color',[.7 .7 .7])
        
        % plot by color
        hold on
        colors = ['rgb']';
        for col = 1:length(average_y)
            plot(average_x{col},average_y{col},['-' colors(col)])
        end
        
        title(datarun.cell_types{cell_type}.name)
        
        drawnow
        
    end

end


if nargout>=1
    plot_fig = params.figure;
end

