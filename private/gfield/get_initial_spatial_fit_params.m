function spatial_params = get_initial_spatial_params(sta, varargin)
% spatial_params = get_initial_spatial_fit_params(sta, varargin)
%   
% provide sta and get initial guess for spatial parameters of Gaussian fit

%% define inputs and parse
p = inputParser;
p.addRequired('sta', @isnumeric);
p.addParamValue('marks_params', [], @isstruct);
p.addParamValue('robust_percent', 0.75, @isnumeric);

p.parse(sta, varargin{:});
marks_params = p.Results.marks_params;
robust_percent = p.Results.robust_percent;

%% body of function

% get the significant stixels to figure out initial conditions of the fit
if isempty(marks_params)
    sig_stixels = significant_stixels(sta);
else
    sig_stixels = significant_stixels(sta, marks_params);
end
% get matrix subscripts to these pixels
[matrix_subscript_i, matrix_subscript_j] = find(sig_stixels);
matrix_subscripts = [matrix_subscript_i, matrix_subscript_j];

% compute the initial center of spatial fit
if size(matrix_subscripts,1) > 1
    initial_center = centroid(matrix_subscripts);
else
    initial_center = matrix_subscripts;
end
    
% comput the intial sd of spatial fit
if size(matrix_subscripts,1) > 1;
    distance_list = ipdm(matrix_subscripts);
    distance_list = triu(distance_list);
    [~, ~, distance_list] = find(distance_list);
    initial_sd = robust_mean(distance_list, robust_percent);
else
    initial_sd = 1;
end


initial_angle = 0;

spatial_params.initial_center = initial_center;
spatial_params.initial_sd = initial_sd;
spatial_params.initial_angle = initial_angle;


