function stacks = stack_align(stacks, input_index, base_index, varargin)
% STACK_ALIGN    Align two stack coordinates to each other using CPSELECT
% usage: stacks = stack_align(stacks, input_index, base_index, opts)
%
% inputs:   stacks          Cell array of image stack structs, as described in Proposal.rtf
%           input_index     Index of starting coordinate space
%           base_index      Index of goal/base coordinate space
%
% opts:     input_slice         1               Slice to use from input
%           base_slice          1               Slice to use from goal/base
%           make_tform          false           Whether to make transform once points are selected
%           tform_type          'polynomial'    Transform type, see also maketform
%           make_array_tform    false           Whether to generate transform to array coordinates,
%                                                   assuming the new transform provides a route and 
%                                                   no transform to array coordinates yet exists.
%           array_tform_type    'polynomial'
%           verbose             true
%
% outputs: Saves control points to stack struct, as well as any newly
% generated transforms.  Also regenerates transform routes.
%
% See also: build_tform_routes build_tform_from_routes array_align
%
% 2010-08 phli
%

opts = inputParser;
opts.addParamValue('input_slice', 1);
opts.addParamValue('base_slice',  1);
opts.addParamValue('make_tform', false);
opts.addParamValue('tform_type', 'polynomial');
opts.addParamValue('make_array_tform', false);
opts.addParamValue('array_tform_type', 'polynomial');
opts.addParamValue('verbose', true);
opts.parse(varargin{:});
opts = opts.Results;

input_index = parse_index(input_index);
base_index  = parse_index(base_index);

istack = get_stack(stacks, input_index);
bstack = get_stack(stacks, base_index);

ips = get_input_points(istack, base_index);
bps = get_base_points( istack, base_index);

islice = get_slice(istack, opts.input_slice);
bslice = get_slice(bstack, opts.base_slice);

% Cut alpha channels if necessary
if size(islice,3) == 4, islice = islice(:,:,1:3); end
if size(bslice,3) == 4, bslice = bslice(:,:,1:3); end

% Unfortunately, CPSELECT will not take empty IN_CPS/BASE_CPS so we have to
% work around with some control flow...
if ~isempty(ips) && ~isempty(bps)
    [ips, bps] = cpselect(islice, bslice, ips, bps, 'Wait', true);
else
    [ips, bps] = cpselect(islice, bslice, 'Wait', true);
end


% Save
istack.input_points{base_index{:}} = ips;
istack.base_points{base_index{:}}  = bps;
stacks{input_index{:}} = istack;


% There are various ways these extra finishing steps can generate errors.
% If they generate errors, we don't want to lose all the work the user did
% in selecting control points.  So we catch errors and convert them into
% mere warnings.
stacks = try_warn(@stack_align_extras, {stacks, input_index, base_index, opts}, 'all', 'Error in generating transforms', stacks);



function stacks = stack_align_extras(stacks, input_index, base_index, opts)
istack = stacks{input_index{:}};

% By default, also make the transform for these control points
if opts.make_tform && isempty(get_stack_tforms(istack, base_index))
    stacks = make_stack_tform(stacks, input_index, base_index, opts.tform_type);
end


% By default, also build the transform to the array if this is now possible
if opts.make_array_tform && isempty(get_stack_tforms(istack, 'array')) && ~isempty(get_tform_routes(istack, 'array'))
    if opts.verbose, disp('Building transform to array coordinates...'); end
    stacks = build_tform_from_routes(stacks, input_index, 'array');
end