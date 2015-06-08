% AKHeitman 2014-04-27
% Concise way to get ROI coord given center and length


function [ROIcoord] = ROI_coord(ROI_length, center, stimsize)
% ROI_length preferably odd
% Stimsize.width (corresponds to x)
% Stimsizse.height  (corresponds to y)
% center.x_coord, center.y_coord


% ROI size cannot exceed screen dimensions
if ROI_length > stimsize.width || ROI_length > stimsize.height
    ROI_length = min(stimsize.width, stimsize.height);
end

% Initial coord
modvec = -floor(ROI_length/2) : floor(ROI_length/2) ;
xdim = center.x_coord + modvec;
ydim = center.y_coord + modvec;

% Correct if ROI falls off the screen
if min(xdim)<1
	xdim = 1:ROI_length;
end
if min(ydim)<1
	ydim = 1:ROI_length;
end
if max(xdim) > stimsize.width
	xdim = (stimsize.width - ROI_length + 1):stimsize.width ;
end
if max(ydim) > stimsize.height
	ydim = (stimsize.height - ROI_length + 1):stimsize.height ;
end

% Assign the final answer
ROIcoord.xvals =  xdim;
ROIcoord.yvals  = ydim;


end
