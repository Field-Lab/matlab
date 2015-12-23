function returned_sta = sta_color_transformed(sta, transform, varargin)
% sta_color_transformed     Applies transform to color dimension of stas
%
% usage:  datarun = sta_color_transformed(sta, transform, params)
%
% arguments:      sta - a single sta (4D array)
%           transform - a 3x3 matrix 
%              params - struct of optional parameters (see below)
%
% outputs:     
%        returned_sta - color transformed STA
%
%
% optional fields in params, their default values, and what they specify:
%   currently params is just a placeholder for future params to be added
% 
%
% GDF 2008-10-1


% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addRequired('sta', @isnumeric);
p.addRequired('transform', @isnumeric);
p.addParamValue('verbose', false,@islogical)

% parse inputs
p.parse(sta, transform, varargin{:});
verbose = p.Results.verbose;

if verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end

% error handling
if length(size(sta)) ~= 4
    error('sta is not a 4D array')
end

if size(sta,3) == 1
    warning('sta is black and white')
end


% APPLY TRANSFORM TO STA

% reshape STA so that matrix multiplication can be used to transform the
% colors, then reshape back into an sta.
temp_sta = reshape(permute(sta, [1 2 4 3]), [], 3);
temp_sta = temp_sta * transform;
temp_sta = permute(reshape(temp_sta, size(sta,1),size(sta,2),size(sta,4),size(sta,3)), [1 2 4 3]);

returned_sta = temp_sta;

% display how long it took
if verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end
