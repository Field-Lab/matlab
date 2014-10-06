function varargout = electrodes_in_stack(datarun, varargin)
% ELECTRODES_IN_STACK    Determine which electrodes are within stack image bounds
% usage: [electrodes, electrodes_tf, edges_xy] = electrodes_in_stack(datarun, stack_index)
% usage: [electrodes, electrodes_tf, edges_x, edges_y] = electrodes_in_stack(datarun, stack_index)
%
% outputs:  electrodes      The indices of the electrodes within the stack image bounds
%           electrodes_tf   Boolean true/false for each electrode, indicating whether inside bounds
%           edges_xy        Nx2 coordinates giving polygon of image bounds
%           [edges_x, edges_y]  Same as edges_xy, given in two separate vectors for convenience
%
% 2010-08 phli
%

stack_index = parse_stack_index(varargin{:});
stack = get_stack(datarun, stack_index);
positions = datarun.ei.position;
[elecs, in, edgesxt, edgesyt] = points_in_stack(stack, positions);

if nargout == 0
    plot_electrodes(positions, {elecs}, 'scale', 0.5);
    hold on;
    plot(edgesxt, edgesyt);
else
    varargout{1} = elecs;
    varargout{2} = in;
end

if nargout == 3
    varargout{3} = [edgesxt' edgesyt'];
end

if nargout > 3
    varargout{3} = edgesxt;
    varargout{4} = edgesyt;
end