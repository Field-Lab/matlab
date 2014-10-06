function roi = ROI_from_contours(sta_conts, params)
% ROI_FROM_CONTOURS     Generate a ROI (region of interest) from STA contours
%
% usage:  datarun = my_function(datarun, arg1, params)
%
% arguments:     sta_conts - contours in standard format
%                   params - struct of optional parameters (see below)
%
% outputs:     roi - logical matrix containing the ROI
%
%
% optional fields in params, their default values, and what they specify:
%
% roi_size          [320 640]   size of ROI matrix
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.roi_size = [320 640];

% combine user and default parameters
params = default_params( defaults, params);






% intilialize variable
roi = zeros(params.roi_size);

% add in each contour
for bb = 1:length(sta_conts)
    
    if sta_conts(bb).hole
        cnt_sign = -1;
    else
        cnt_sign = 1;
    end
    
    roi = roi + cnt_sign * roipoly(roi, sta_conts(bb).x, sta_conts(bb).y);
end

