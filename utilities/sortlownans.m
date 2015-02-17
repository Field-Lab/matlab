function [S, inds] = sortlownans(A)
% SORTLOWNANS   Runs sort but with NaNs sorting to the low end instead of the high end
%
% 2012-09 phli
%

B = A;
B(isnan(B)) = -Inf;
[~, inds] = sort(B);
S = A(inds);