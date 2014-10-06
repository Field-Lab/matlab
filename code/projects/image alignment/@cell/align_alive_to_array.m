function stacks = align_alive_to_array(stacks, alive_index, datarun_or_array_id, varargin)
% ALIGN_ALIVE_TO_ARRAY    Quick alignment method by clicking corners
% usage: stacks = align_alive_to_array(stacks, alive_index, datarun_or_array_id, opts)
%
% inputs:   stacks                  A cell array of stack structs, as described in Proposal.rtf
%           alive_index             The index of the stack, usually {2,1} or {2,2}
%           datarun_or_array_id     Either a datarun struct or simply the array id, e.g. 504
%
% opts:     slice   1               Which slice from the stack to use for aligning
%           For other opts and more info see also align_alive_to_array
%
% outputs: Saves into the stack struct the control points and transform to
% array and array_image coordinates
%
% 2010-08 phli
%

opts = inputParser;
opts.addParamValue('slice', 1);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

alive_index = parse_stack_index(alive_index);
stack = get_stack(stacks, alive_index);
im = get_slice(stack, opts.slice);
[tform_array, array_corners, ai] = align_alive_to_array(im, datarun_or_array_id, unmatched);


% Save
corner_electrode_positions = ai.positions(ai.corner_electrodes,:);
stack.input_points{1} = array_corners;
stack.base_points{1} = corner_electrode_positions;
stack.tforms{1} = tform_array;
stack.tforms_inv{1} = fliptform(tform_array);

% Now make control points for transform to array_image as well
array_image_corners = tformfwd(ai.T_array_to_array_image, corner_electrode_positions);
stack.input_points{1,2} = array_corners;
stack.base_points{1,2} = array_image_corners;

% Save
stacks{alive_index{:}} = stack;

% Fill in transforms to/from array_image based on the control points added above
stacks = make_stack_tform(stacks, alive_index, 'array_image', 'projective');