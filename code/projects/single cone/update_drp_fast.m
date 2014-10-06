%
% helper function for generate_clumped_mosaic_using_drp
%
% the drp being updated is computed from targetColor to referenceColor
% since exactly two labels have changed in c, we simply need to update the
% distance matrix from the last iteration (allDistances) appropriately
% and then recompute the density recovery profile
function [drp,counts] = update_drp_fast(c,allDistances,numBins,binWidth,targetColor,referenceColor)

% compute bin centers with an extra bin at the end
edges = (0:1:numBins)*binWidth;
edges = [-Inf edges Inf];

% get area of each annulus
areas = pi*binWidth^2*(2*(1:numBins)-1);

% choose appropriate distances
if targetColor == referenceColor
    dists = squareform(allDistances(c{referenceColor},c{targetColor}));
    areas = areas * 0.5;
else
    dists = reshape(allDistances(c{referenceColor},c{targetColor}),[],1)';
end

% bin stuff up
counts = histc(dists(dists<Inf),edges);

% kill first and last bins
counts = counts(2:end-2);

% divide by the area and the normalization for the number of points
% also scale by count(target)/count(reference)
drp = (counts./(areas*length(c{targetColor})))*(length(c{targetColor})/length(c{referenceColor}));