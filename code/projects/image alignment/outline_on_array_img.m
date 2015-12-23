function outline_on_array_img(stacks, index, islice)
% OUTLINE_ON_ARRAY_IMG    Show outline of stack image bounds on array image
% usage: outline_on_array_img(stacks, index, islice)
%
% inputs:   stacks  Cell array of image stack structs as described in Proposal.rtf
%           index   Index of stack to outline
%           islice  Slice to use for outline, although this shouldn't
%                       matter as slices within our stacks should have
%                       uniform dimensions.  Defaults to 1.
%
% 2010-08 phli
%

if nargin < 3
    islice = 1;
end


aii = parse_stack_index('array_image');

index = parse_stack_index(index);
stack = get_stack(stacks, index);
array_image_tform = stack.tforms{aii{:}};
if isempty(array_image_tform)
    error('No transform to array image!');
end

[edgesx, edgesy] = stack_edges(stack, islice);
[edgesxt, edgesyt] = tformfwd(array_image_tform, edgesx, edgesy);

imshow(get_slice(stacks{aii{:}}));
hold on;
plot(edgesxt, edgesyt);