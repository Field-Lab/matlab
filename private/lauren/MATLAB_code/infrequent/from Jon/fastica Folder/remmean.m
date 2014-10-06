function [newVectors, meanValue] = remmean (vectors);

% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
%
% This function is needed by FASTICA and FASTICAG

% 15.3.1998

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));