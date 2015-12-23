function [square_equivalent,pixel_height,pixel_width] = compute_pixel_size(datarun, varargin)
% compute_pixel_size     compute the size of pixels in physical coordinates
%
%   uses the tranformation from monitor space to array space.  because the monitor projection is not perfectly
%   uniform, pixels at different locations have different sizes when projected onto the array.  in this code, 
%   the default behavior is to use a pixel at the center of the monitor to compute the size.
%
%
% usage:  [square_equivalent,pixel_height,pixel_width] = compute_pixel_size(datarun, <params>)
%
% arguments:     datarun - requires that datarun.piece.T_monitor_to_array exists
%                           see compute_monitor_to_array_transformation
%               <params> - struct or list of optional parameters (see below)
%
% outputs:     square_equivalent - sidelength of a square with area equal to one pixel (microns)
%                   pixel_height - height of one pixel (microns)
%                    pixel_width - width of one pixel (microns)
%
%
% optional params, their default values, and what they specify:
%
% monitor_location       	[datarun.stimulus.monitor_x datarun.stimulus.monitor_y]/2
%                               what point in monitor space to use to compute the pixel size
%                               defaults to the center of the monitor
%
% See also compute_monitor_to_array_transformation
%
% 2010-06  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('monitor_location', [datarun.stimulus.monitor_x datarun.stimulus.monitor_y]/2);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% get the transformation
T = datarun.piece.T_monitor_to_array;

% create the coordiantes for a horizontal segment one pixel long, and a vertical segment one pixel high
% by default, they are located at the center of the stimulus monitor
monitor_segment_horizontal = [ params.monitor_location(1) + [0 1]'  params.monitor_location(2) + [0 0]'];
monitor_segment_vertical = [ params.monitor_location(1) + [0 0]'  params.monitor_location(2) + [0 1]'];

% transform these segments into array space, where the units are microns
array_segment_horizontal = tformfwd(T,monitor_segment_horizontal);
array_segment_vertical = tformfwd(T,monitor_segment_vertical);

% compute the length of the segments in array space
pixel_width = norm(array_segment_horizontal(1,:)-array_segment_horizontal(2,:));
pixel_height = norm(array_segment_vertical(1,:)-array_segment_vertical(2,:));

% compute the geometric mean
square_equivalent = sqrt(pixel_width*pixel_height);


