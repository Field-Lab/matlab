function indexed = index_masks(masks, indexes)
% INDEX_MASKS   Convert masks to indexed pixel bitmap as specified
% usage: indexed = index_masks(masks, indexes)
%
% inputs: MASKS     Cell array where each element is an MxN binary mask
%                   representing a single cone receptive field.
%         INDEXES   Cell array with series of lists of masks to combine
%                   into given index.  E.g. to put a set of masks combined
%                   into index 4 and a set into index 5:
%                       indexes = {};
%                       indexes{4} = 1:100;
%                       indexes{5} = 101:110;
%
% output: INDEXED   HEIGHT x WIDTH pixel bitmap with integer values
%                   indicating the indices.  Should work directly with
%                   IMAGE() (treats it as an indexed graphic like a GIF or 
%                   PNG), unless indices run to greater than 255, in
%                   which case you could try IMAGESC();
%
% See also: MAKE_VORONOI_MASKS, IMAGE, IMAGESC
% 
% phli 2011-01
%

if nargin < 2
    indexes = num2cell(1:length(masks));
end

notempty = detect(masks, @(e) (~isempty(e)));

indexed = zeros(size(notempty, 1), size(notempty, 2));
for i = 1:numel(indexes)
    imasks = masks(indexes{i});
    mask = combine_masks(imasks);
    indexed(mask) = i;
end