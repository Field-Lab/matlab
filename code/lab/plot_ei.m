function elec_handles = plot_ei(datarun, cell_id, varargin)
% plot_ei     plot an electrophysiological image
%
% usage:  plot_ei(datarun, cell_id, varargin)
%
% arguments:     datarun - datarun struct
%                cell_id - cell id of the EI to plot
%               varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% foa               gca            	figure or axes to plot in. if 0, make new figure.
%                                      	if empty, don't plot.  if -1, plot in current.
% frame_number      0               frame number to plot.  if 0, plot each electrode's maximal frame%
% coordinates       'vision ei'     which coordinate frame to plot in
%                                       'vision ei'
%                                       'array'
%                                       'array image'
%                                       'monitor'
%                                       'sta'
%
% all other parameters are passed to the low level function plot_ei_.  see that function for details.
%
%
%
% 2009      greschner
% 2010-01   gauthier, added option for coordinates
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% some provided arguments might be for the low level function and should be ignored.
p.KeepUnmatched = true;

% specify list of optional parameters
p.addParamValue('foa', gca);
p.addParamValue('frame_number', 0);
p.addParamValue('coordinates', 'vision ei',@(x)any(strcmpi(x,{'vision ei','array','array image','monitor','sta'})));
p.addParamValue('sta_udflip', 0);
p.addParamValue('sta_lrflip', 0);
p.addParamValue('stack', {});
% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


if ~isempty(params.stack)
    switch class(params.stack)
        case 'cell'
            stack = datarun.stacks{params.stack{:}};
        case 'struct'
            stack = params.stack;
    end
    outline = stack_edges(stack);
else
    outline = [];
end


% transform electrode positions...
switch params.coordinates
    case 'vision ei'  % to vision EI coordinates
        array_info = load_array_info(datarun,2);
        positions = tformfwd(array_info.T_array_to_vision_ei,datarun.ei.position);
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
            outline = tformfwd(array_info.T_array_to_vision_ei, outline);
        end
    case 'array'  % to array coordinates
        positions = datarun.ei.position;
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
        end
        
    case 'array image'  % to array image coordinates
        array_info = load_array_info(datarun,2);
        positions = tformfwd(array_info.T_array_to_array_image,datarun.ei.position);
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1,2}, outline);
        end
        
    case 'monitor' % to monitor coordinates
        T =  fliptform(datarun.piece.T_monitor_to_array);
        positions = tformfwd(T,datarun.ei.position);
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
            outline = tformfwd(T, outline);
        end
        
    case 'sta' % to STA coordinates
        T_array_to_monitor = fliptform(datarun.piece.T_monitor_to_array);
        T_monitor_to_STA = fliptform(coordinate_transform(datarun, 'monitor'));
        T = maketform('composite',[T_monitor_to_STA T_array_to_monitor]);
        positions = tformfwd(T, datarun.ei.position);
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
            outline = tformfwd(T, outline);
        end
        if params.sta_udflip
            positions(:,2)=datarun.stimulus.field_height-positions(:,2);
            if ~isempty(outline)
                outline(:,2) = datarun.stimulus.field_height - outline(:,2);
            end
        end
        if params.sta_lrflip
            positions(:,1)=datarun.stimulus.field_width-positions(:,1)+.5;
            if ~isempty(outline)
                outline(:,1) = datarun.stimulus.field_width - outline(:,1);
            end
        end
    otherwise
        error('coordinates spec ''%s'' not recognized.',params.coordinates)
end
 

% call low level function
elec_handles = plot_ei_(get_ei(datarun,cell_id), positions, params.frame_number, p.Unmatched, 'stack_outline', outline);
if nargout < 1; clear elec_handles; end