function stacks = build_tform_from_routes(stacks, start_index, goal_index, varargin)
% BUILD_TFORM_FROM_ROUTES    Create transform based on determined routing
% usage: stacks = build_tform_from_routes(stacks, start_index, goal_index, opts)
%
% inputs:   stacks       Cell array of stack structs as described in Proposal.rtf
%           start_index  Stack index of starting coordinates, e.g. {2,1} by convention for alive coarse montage
%           goal_index   Stack index of ending coordinates, e.g. {1,1} by convention for array coordinates
%
% opts:     method  'composite'     Method to use to contruct transform.
%                                       Currently only 'composite' is accepted.
%           routes  []              Explicitly selected routes; usually not
%                                       needed
%
% outputs: saves new transforms in stacks{start_index{:}}.tforms{goal_index{:}} and .tforms_inv
%
% See also: 
%
% 2010-08 phli
%

opts = inputParser;
opts.addParamValue('method', 'composite');
opts.addParamValue('routes', []);
opts.parse(varargin{:});
opts = opts.Results;

start_index = parse_stack_index(start_index);
goal_index  = parse_stack_index(goal_index);
stack = get_stack(stacks, start_index);

% Check that the route actually exists
if isempty(get_tform_routes(stack, goal_index))
    error('No appropriate route!');
end

% Initialize parameters for while loop
current_stack = stack;
current_index = cell2mat(start_index);
step_number = 1;

% Collect route to array
while ~all(current_index == cell2mat(goal_index))
    current_routes = get_tform_routes(current_stack, goal_index);
    if length(opts.routes) >= step_number && ~isempty(opts.routes(step_number))
        current_route = current_routes{opts.routes(step_number)};
    else
        current_route = current_routes{1};
    end

    % Add forward and inverse transforms to the chains
    next_index = num2cell(current_route);
    [current_tform, current_tform_inv] = get_stack_tforms(current_stack, next_index);
    if step_number == 1
        tforms     = current_tform;
        tforms_inv = current_tform_inv;
    else
        tforms(end+1)     = current_tform; %#ok<AGROW>
        tforms_inv(end+1) = current_tform_inv; %#ok<AGROW>
    end
    
    % Next loop
    current_stack = get_stack(stacks, next_index);
    current_index = current_route;
    step_number = step_number + 1;
end

% Convert chain of transforms into a single equivalent
switch opts.method
    case 'points'
        error('Not implemented yet!');

    case 'composite'
        goal_tform     = maketform('composite', fliplr(tforms));
        goal_tform_inv = maketform('composite', tforms_inv);
        
    otherwise
        error('Unrecognized transform method');
end

% Save
stack.tforms{goal_index{:}} = goal_tform;
stack.tforms_inv{goal_index{:}} = goal_tform_inv;
stacks{start_index{:}} = stack;