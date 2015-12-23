function routes = get_tform_routes(stack, varargin)
% GET_TFORM_ROUTES
% usage: routes = get_tform_routes(stack, index)
%
% 2010-08 phli
%

index = parse_stack_index(varargin{:});

if isfield(stack, 'tform_routes') && all(size(stack.tform_routes) >= [index{:}])
    routes = stack.tform_routes{index{:}};
else
    routes = [];
end