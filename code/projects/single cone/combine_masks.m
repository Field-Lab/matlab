function combined = combine_masks(masks)
% COMBINE_MASKS     Combine pixel masks
% 
% phli 2011-01
%

nonemptymask = detect(masks, @(m) (~isempty(m)));
combined = false(size(nonemptymask));

for i = 1:numel(masks)
    mask = masks{i};
    if isempty(mask)
        warning(['Empty mask: ' num2str(i)]);
        continue;
    end
    
    combined = combined | mask;
end