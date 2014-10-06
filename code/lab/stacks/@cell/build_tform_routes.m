function stacks = build_tform_routes(stacks, varargin)
% BUILD_TFORM_ROUTES    Build chains of transforms to reach given coordinates
% usage: stacks = build_tform_routes(stacks, goal)
%
% inputs: stacks    Cell array of stack structs as described in Proposal.rtf
%         goal      The index for the coordinates to try to reach.  Most
%                       commonly we are interested in generating transform
%                       chains that end in array coordinates, i.e. index {1,1}.
%
% outputs: saves routing information into the stack structs, in the field tform_routes
%
% 2010-08 phli
%

% If no goal is given, recursively run for every possibility!
if nargin < 2
    for i = 1:numel(stacks)
        [x y] = ind2sub(size(stacks), i);
        stacks = build_tform_routes(stacks, {x,y});
    end
    
    return;
end

goal = parse_stack_index(varargin{:});

% Get 1-D indices of elements of stacks that aren't empty
toroute = cell2mat(collect(stacks(:), @(elem) ~isempty(elem)))';

% But we don't need to route the goal to itself; that's already done
igoal = sub2ind(size(stacks), goal{:});
toroute(igoal) = false;
last_routed = igoal;

% Start from goal node and iteratively add routes
while any(toroute)
    routed = [];

    % Iterate through the stacks that still need to be routed
    for i = find(toroute)
        istack = get_stack(stacks, i);
        istack.tform_routes{goal{:}} = {};
        
        % Check whether there is any transform from this stack to a stack
        % that has already been routed
        if ~isfield(istack, 'tforms')
            continue;
        end
        for j = last_routed
            % Convert the 1-D index j back into a 2-D subscript
            [x y] = ind2sub(size(stacks), j);
            
            if all(size(istack.tforms) >= [x y]) && ~isempty(istack.tforms{x,y})
                % There is a transform from here to an already routed node,
                % so save this transform as the routing for this stack.
                %
                % More concise method istack.tform_routes{goal{:}}{end+1} = ... failed; seems to be a MatLab subscripting bug
                routes = istack.tform_routes{goal{:}};
                routes{end+1} = [x y]; %#ok<AGROW>
                istack.tform_routes{goal{:}} = routes;
                stacks{i} = istack;

                % Note that this one has been routed on this iteration
                routed(end+1) = i; %#ok<AGROW>
            end
        end
    end
    
    % We didn't add anything this round, so time to give up
    if isempty(routed)
        break;
    else
        % Mark off the ones that no longer need to be routed
        toroute(routed) = false;
        
        % Prep for next iteration
        last_routed = unique(routed);
    end
end