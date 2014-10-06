function [slice,sample_points] = matrix_slice(matrix,segment,varargin)
% matrix_slice     Identify values in a matrix along a specified segment.
%
% note: any points which fall outside the bounds of the matrix will be NaN
%
% usage:  slice = matrix_slice(matrix,segment,varargin)
%
% arguments:   matrix - XxYxC matrix.  values are sampled in the xy plane along the specified segment.
%             segment - coordinates for the segment endpoints, [x_start y_start; x_end y_end]
%            varargin - struct or list of optional parameters (see below)
%
% outputs:      slice - NxC matrix of values in the original matrix
%       sample_points - the points that were actually sampled
%
%
% optional params, their default values, and what they specify:
%
% samples           100         how many samples to take along the specified segment
% extension         0           how much to extend the segment beyond the specified endpoints
%                               units are the length of the segment
%                               e.g. extension = 0.6 means the total slice will
%                                  	be (1 + 0.6 + 0.6) times as long as the specified segment
%
%
% 2009-04 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('samples', 100, @(x)(x>0 && round(x)==x)); % must be postive integer
p.addParamValue('extension', 0, @(x)(x>=0))

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% verify matrix is not empty
if isempty(matrix)
    slice = [];
    sample_points = [];
    return
end


% GET COORDINATES FOR EACH SAMPLE POINT

% rename
pt_1 = segment(1,1:2);
pt_2 = segment(2,1:2);

% extrapolate points to include extension factor
y_range = pt_2(2)-pt_1(2);
x_range = pt_2(1)-pt_1(1);
pt_1 = pt_1 - params.extension * [x_range y_range];
pt_2 = pt_2 + params.extension * [x_range y_range];

% generate coordinates
x_range = pt_2(1)-pt_1(1);
y_range = pt_2(2)-pt_1(2);


% interpolate points between the coordinates

if x_range == 0
    % handle vertical segment
    x_points = repmat(pt_1(1),1,params.samples);
else
    x_points = [pt_1(1):x_range/(params.samples-1):pt_2(1)];
end

if y_range == 0
    % handle horizontal segment
    y_points = repmat(pt_1(2),1,params.samples);
else
    y_points = [pt_1(2):y_range/(params.samples-1):pt_2(2)];
end





% TAKE THE SAMPLES

% initialize the slice variable with NaNs
slice = nan(params.samples, size(matrix,3));

% identify which points are in the bounds of the matrix
in_range = intersect(find(  round(x_points)>=1 & round(x_points)<=size(matrix,2)  ),...
    find(  round(y_points)>=1 & round(y_points)<=size(matrix,1)  ));

% collect matrix values from points which are in range
for nn = in_range
    slice(nn,:) = matrix(round(y_points(nn)),round(x_points(nn)),:);
end




% put points into a single set of coordinates
sample_points = [x_points' y_points'];
