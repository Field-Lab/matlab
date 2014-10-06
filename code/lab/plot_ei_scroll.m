function varargout = plot_ei_scroll(datarun,cell_id, varargin)
% plot_ei_scroll     plot an EI with a scroll bar
%
% usage:  fig = plot_ei_scroll(datarun,cell_id, <params>)
% usage: [fig, output_data] = plot_ei_scroll(datarun, cell_id, 'output_data', appdataname)
%
% arguments:     datarun - datarun struct
%                cell_id - which cell id
%               <params> - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
%
% coordinates       'vision ei'     which coordinate frame to plot in
%                                       'vision ei'
%                                       'array'
%                                       'array image'
%
%
% all other parameters are passed to plot_ei_scroll_ and any parameters
% not accepted there are passed to plot_ei_.  see those functions for options.
%
%
% examples:
%  
%      plot_ei_scroll(datarun,cell_id)
%      plot_ei_scroll(datarun,cell_id,'figure',10)
%      plot_ei_scroll(datarun,cell_id,'start_frame',1)
%      plot_ei_scroll(datarun,cell_id,'neg_color','k')
%
%
%
% 2010-03  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% some provided arguments might be for the low level function and should be ignored.
p.KeepUnmatched = true;

% specify list of optional parameters
p.addParamValue('coordinates', 'vision ei', @(x) any(strcmpi(x,{'vision ei','array','array image','stack'})));
p.addParamValue('axon', []);
p.addParamValue('modify_ei', []);
p.addParamValue('stack', {});
p.addParamValue('output_data', []);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


outline = [];
if ~isempty(params.stack)
    switch class(params.stack)
        case 'cell'
            stack = datarun.stacks{params.stack{:}};
        case 'struct'
            stack = params.stack;
    end

    outline = stack_edges(stack);
end


% get axon from params
axon = params.axon;


% transform electrode positions...
positions = datarun.ei.position;
switch params.coordinates
    case 'vision ei'  % to vision EI coordinates
        array_info = load_array_info(datarun,2);
        positions = tformfwd(array_info.T_array_to_vision_ei,positions);
        if ~isempty(axon)
            axon = tformfwd(array_info.T_array_to_vision_ei,axon);
        end
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
            outline = tformfwd(array_info.T_array_to_vision_ei, outline);
        end
        
    case 'array'  % to array coordinates
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1}, outline);
        end
        
    case 'array image'  % to array image coordinates
        array_info = load_array_info(datarun,2);
        positions = tformfwd(array_info.T_array_to_array_image,positions);
        if ~isempty(axon)
            axon = tformfwd(array_info.T_array_to_array_image,axon);
        end
        if ~isempty(outline)
            outline = tformfwd(stack.tforms{1,2}, outline);
        end
        
    case 'stack'
        ai = parse_stack_index('array');
        array_tform = stack.tforms_inv{ai{:}};
        positions = tformfwd(array_tform, positions);
        if ~isempty(axon)
            axon = tformfwd(array_tform, axon);
        end
        
    otherwise
        error('coordinates spec ''%s'' not recognized.',params.coordinates)
end


% get the ei
ei = get_ei(datarun,cell_id);

% Run modification code
if ~isempty(params.modify_ei) && isa(params.modify_ei, 'function_handle')
    ei = params.modify_ei(ei, datarun, cell_id);
end

% note the frame of the spike
spike_frame = datarun.ei.nlPoints;

% set up title text
prefix = sprintf('cell id %d',cell_id);

% append cell type, if is exists
ct = find_cell_types(datarun,cell_id);
if ct > 0
    prefix = [prefix sprintf(' (%s)',datarun.cell_types{ct(1)}.name)];
end

if isfield(datarun, 'sampling_rate')
    sampling_rate = datarun.sampling_rate;
else
    sampling_rate = [];
end

% call low level function
fig = plot_ei_scroll_(ei, positions,'spike_frame',spike_frame,'prefix',prefix,'axon',axon, 'sampling_rate', sampling_rate, 'stack_outline', outline, p.Unmatched);
if strcmp(params.coordinates, 'stack')
    set(gca, 'YDir', 'reverse');
    axis image;
end

if nargout > 0
    varargout{1} = fig;
end

% Some tricks to get data out of the figure as it closes
if ~isempty(params.output_data) && nargout > 1
    dataout = waitfordata(fig, params.output_data);
    varargout = [varargout dataout];
end