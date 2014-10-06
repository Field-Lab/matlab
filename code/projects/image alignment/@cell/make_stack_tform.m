function stacks = make_stack_tform(stacks, input_index, base_index, tform_type)
% MAKE_STACK_TFORM    Build transform from existing control points
% usage: stacks = make_stack_tform(stacks, input_index, base_index, tform_type)
%
% inputs:   stacks          Cell array of image stack structs as described in Proposal.rtf
%           input_index     Starting coordinates
%           base_index      Goal coordinates
%           tform_type      Type of tform; see also maketform
%
% Will look in stack struct for control points, as returned from
% stack_align or array_align.  Uses control points to build transform and
% save it to stacks.
%
% See also: MAKETFORM
%
% 2010-08 phli
%

if nargin < 4
    tform_type = 'polynomial';
end

input_index = parse_stack_index(input_index);
base_index  = parse_stack_index(base_index);
istack = get_stack(stacks, input_index);

ips = get_input_points(istack, base_index);
bps = get_base_points( istack, base_index);
tform = cp2tform(ips, bps, tform_type);

% Different transform types result in differently filled in transform structs...
if ~isempty(tform.forward_fcn)
    istack.tforms{base_index{:}}     = tform;
    istack.tforms_inv{base_index{:}} = fliptform(tform);
else
    tform_inv = cp2tform(bps, ips, tform_type);
    istack.tforms{base_index{:}}     = fliptform(tform_inv);
    istack.tforms_inv{base_index{:}} = fliptform(tform);
end

% Save
stacks{input_index{:}} = istack;

% Regenerate transform routes
stacks = build_tform_routes(stacks);