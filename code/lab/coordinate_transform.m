function tform = coordinate_transform(datarun, output, varargin)
% coordinate_transform     create transformation between stimulus coordinate spaces, e.g. stixels to pixels
%
%   input space is STA coordinates, either normal or scaled up 
%   output space can be monitor coordinates, STA coordinates, or "rect" (a user defined rotated rectangle, useful
%       for plotting figures that are zoomed in on an interseting collection of RFs)
%
%
% usage:  tform = coordinate_transform(datarun, output, varargin)
%
% arguments:  datarun - datarun struct with stimulus field
%              output - what space to transform to
%                           'monitor' - 
%                           'sta' - 
%                           'rect' - 
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% input          	'sta'    	input space
%                                   'sta' - stixels in the STA
%                                   'sta scaled' - scaled up STA
%
%
%  if input == 'sta scaled'
%       scale      	1         	how much the STA is scaled up
%
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('input', 'sta',@(x) any(strcmpi(x,{'sta','sta scaled'})));
p.addParamValue('scale', 1,    @(x) mod(x,0)==x);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% make transformation
% this transformation will convert from STA coordinates to the desired plotting coordinates


% the function gets numerical coordinates for the upper left and lower right corners of the STA
% the actual numbers depend on what the input and output spaces are


% get points from the input
switch params.input
    case 'sta'
        % input is in STA coordinates
        input_points = [1 1; datarun.stimulus.field_width datarun.stimulus.field_height];
    case 'sta scaled'
        nx = params.scale;
        ny = params.scale;
        input_points = [(0.5 + nx/2) (0.5 + ny/2); ...
            (datarun.stimulus.field_width*nx - nx/2 + 0.5) (datarun.stimulus.field_height*ny - ny/2 + 0.5)];
end
    
% set output coordinates
switch output
    case 'sta'
        base_points = [1 1; datarun.stimulus.field_width datarun.stimulus.field_height];
        bg_size = [datarun.stimulus.field_width datarun.stimulus.field_height];
    case 'monitor'
        nx = datarun.stimulus.stixel_width;
        ny = datarun.stimulus.stixel_height;
        base_points = [(datarun.stimulus.x_start + nx/2 + 0.5) (datarun.stimulus.y_start + ny/2 + 0.5) ;...
            (datarun.stimulus.x_end - nx/2 + 0.5) (datarun.stimulus.y_end - ny/2 + 0.5)];
        bg_size = [datarun.stimulus.monitor_x datarun.stimulus.monitor_y];
    case 'rect'
        error('not yet implemented')
    otherwise
        error(['Coordinate space not recognized: ' output])
end


% compute transform
switch version('-release')
    case '2007a'
        tform = cp2tform(input_points, base_points, 'linear conformal');
    otherwise
        tform = cp2tform(input_points, base_points, 'nonreflective similarity');
end


