function buf = stavect2visbuf(stavect, width, height)
% Transform Matlab STA vector into Vision STA buffer format vector
% 
% Normal Matlab STA format is Y-X-C-T in separate dimensions.  Vision 
% buffer format is C-X-Y dimensions collapsed in that order in a single 
% vector for each time frame.
%
% A particular case where combining Matlab STA formatted data with Vision
% buffer formatted data is useful is if the Matlab data can be made sparse.
% Matlab sparse matrices cannot be >2D, so this generally involves making
% the Matlab STA into a vector of Y-X-C collapsed for each time frame.  
% E.g. the Wc matrix from cone finding or the STRF matrices used as input
% to generator signal calculations.
%
% STAVECT2VISBUF is a function to convert Matlab YXC sparse vectors to CXY
% Vision buffer vectors.
%
% 2013-06 phli
%

ind = 1:size(stavect,1);
ind = reshape(ind, height, width, []);
ind = permute(ind, [3 2 1]);
ind = reshape(ind, [], 1);
buf = stavect(ind,:);