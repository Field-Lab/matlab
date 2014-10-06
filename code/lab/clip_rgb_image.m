function rgb_out = clip_rgb_image(rgb, thresh)
% clip_rgb_image     ensure no pixels go outside a certain range
%                       any pixel with an R, G, or B value outside the range is divided down 
%                       but the pixel maintins its R-G-B ratio
%
% usage:  rgb_out = clip_rgb_image(rgb, thresh)
%
% arguments:     rgb - YxXxC matrix
%            thresh - 
%
% outputs:      rgb_out - 
%
%
% 2010-05-13  gauthier
%


% get threshold
switch numel(thresh)
    case 1
        pos_thresh = abs(thresh);
        neg_thresh = -abs(thresh);
    case 2
        % must have one positive element and one negative element
        pos_thresh = max(thresh);
        neg_thresh = min(thresh);
    otherwise
        error('threshold must be  two element vector')
end


% reshape
orig_size = size(rgb);
rgb = reshape(rgb,[],size(rgb,3));

% identify pixels outside the range in the positive direction
max_vals = max(rgb,[],2);
over_max = find(max_vals > pos_thresh);

% normalize them
rgb(over_max,:) = rgb(over_max,:) ./ repmat(max_vals(over_max),1,orig_size(3)) .* repmat(pos_thresh,length(over_max),orig_size(3));



% identify pixels outside the range in the negative direction
min_vals = min(rgb,[],2);
under_min = find(min_vals < neg_thresh);

% normalize them
rgb(under_min,:) = rgb(under_min,:) ./ repmat(min_vals(under_min),1,orig_size(3)) .* repmat(neg_thresh,length(under_min),orig_size(3));



rgb_out = reshape(rgb,orig_size);
