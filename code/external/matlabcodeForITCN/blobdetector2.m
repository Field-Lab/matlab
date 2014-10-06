function [i,j] = blobdetector2(im, side, min_dist,option,th,map)
% Modified from the UCSB ICTN library
%
% 2010-06 phli
%

% function p = blobdetector(im, side, min_dist,option,th,map);
%
% Input
%-------------------------------------------------------------
% im = image to count cells
% side: filter size - choose the diameter (or slightly larger) of the blob (in pixels)
% min_dist : minimum distance between peaks
% option : 1 if the peak to detect is dark
%          0 if the peak to detect is bright (e.g. Topro stained image)
% map : mask of layer
% th : threshold for filter output
      %- suggested to choose -0.05 or 0

% Output
%-------------------------------------------------------------
% p: returns the number of peaks detected.
% resultImaage: under current directory as 
% Result(filter size,min_dist)_p_totalarea

% Works only with gray-scale images. If RGB image is input, the image will
% be converted to gray image


if nargin < 4
    option = 0;
end

if nargin < 5
    th = 0;
end

% if map is not specified, assume a whole image as a map
if nargin < 6
    map = ones(size(im,1),size(im,2));
end

if size(im,3) == 1
    img = im;
else
    img = rgb2gray(im);
end

if option == 0
    img = 255 - img;
end



% apply the filter
ac = lapofgau(img, side);
ac(ac < th) = th;
ac = ac - th;

% find local maxima from filter output
scaling = 1;
[i,j] = find_local_max_2D(ac, [], floor(side/4/scaling), inf, min_dist, [], [1 1], map);



if nargout < 2
    figure; imagesc(ac); axis image; axis off
    figure; imagesc(im); axis image; colormap gray; axis off
    hold on; plot(scaling*j, scaling*i, '.r'); hold off
    area = length(find(map == 255 | map == 1));
    title([num2str(length(i)),' blobs (area:',num2str(area),')']);
end