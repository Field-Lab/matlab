function [dist diff] = getDistances(locations,centers)

% function dist = getDistances(locations,centers)
%
% get all pairwise distances between two lists of points. 
% locations is N by k and centers is C by k. k must be the same, but
% N and C can be different. output is a matrix containing the distances,
% with dimension N by C. also returns the difference matrix as a second
% argument 

%%
% get dimensions
N = size(locations,1);
C = size(centers,1);
% replicate locations in third dimension
locsRep = repmat(locations,[1 1 C]);
% replicate centers in first dimension
tmp(1,:,:) = centers';
centsRep = repmat(tmp,[N 1 1]);
% take all center x location differences
diff = (locsRep - centsRep);
%diff = permute(diff,[1 3 2]);
% square and sum x and y component to get a distance
%dist = sum(diff,3);
dist = sqrt(squeeze(sum(diff.^2,2)));
