function [fit_info , sta, sig_stixels] = fit_sta_spatial(sta, varargin)

% fit_sta.m fits a spatial-temporal-chromatic STA with a differences of
% Gaussians (in space), a difference of a cascade of filters (in time), and
% a 1x3 vector of color weights.  It handles BW and RGB STAs.
%
% fit_info = fit_sta(sta, varargin)
%
% INPUTS:
%   sta             a spatial-temporal-chromatic STA
%
% OPTIONAL INPUTS (HANDLED WITH VARARGIN AND INPUTPARSER):
%
%   sig_stixels         []  a sparce matrix of significant stixels, see 
%                           significant_stixels.m 
%   time_course         []  a time course for identifying initial
%                           conditions for temporal and chromatic fits.
%   marks_params        []  parameters for choosing signficant stixels, see
%                           significant_stixels.m
%   tc_params           []  parameters for choosing time course of STA
%                           See time_course_from_sta.m
%
%   x_dim               []  height of STA, computed from data if empty
%   y_dim               []  width of STA, computed from data if empty
%   num_colors          []  number of color channels in STA
%   frame_number        []  number of STA frames
%
%   initial_center_point_x          []  if empty, computed from STA
%   initial_center_point_y          []  ditto
%   initial_center_sd_x             []  ditto
%   initial_center_sd_y             []  ditto
%   initial_center_rotation_angle   []  ditto
%
%   initial_color_weight_a          1
%   initial_color_weight_b          1
%   initial_color_weight_c          1
%   
%   initial_surround_sd_scale       2.0
%   initial_surround_amp_scale      0.15
%
%   initial_scale_one               []      temporal filter 1 amplitude
%   initial_scale_two               []      temporal filter 2 amplitude
%   initial_tau_one                 []      filter 1 time constant
%   initial_tau_two                 []      filter 2 time constant
%   initial_n_one_filters               6       number of filters in each
%   initial_n_two_filters               6       number of filters in each

%                                           cascade
%   fit_center_point_x              true
%   fit_center_point_y              true
%   fit_center_sd_x                 true
%   fit_center_sd_y                 true
%   fit_center_rotation_angle       true
%
%   fit_surround                    false
%   fit_surround_sd_scale           false
%   fit_surround_sd_amp             false
%
%   fit_scale_one                   true
%   fit_scale_two                   true
%   fit_tau_one                     true
%   fit_tau_two                     true
%   fit_n_one_filters                   false
%   fit_n_two_filters                   false
%
%   fit_color_weight_a              true    unless STA is BW
%   fit_color_weight_b              true    unless STA is BW
%   fit_color_weight_c              true    unless STA is BW
%
%   fit_sig_stixels_one             false   only fit significant stixels
%
%   verbose                         false   plot fit results on each
%                                           iteration
%   optim   {'TolFun',0.001,'Display','off',...         See fminsearch doc
%            'MaxIter', 10000, 'MaxFunEvals', 10000}                        
%
%   OUTPUT:
%       fit_info        A structure that contains all the information about
%                       the fit, including initial conditions
%
%
%   FIT_STA finds fits spatial-temporal-chromatic STAs, where STA is
%   [height, width, colors, time_frames] (a 4-D array). Any parameter of
%   the fit can be provided an initial condition, but if one is not
%   provided, then the initial conditions are computed (sometimes ad hoc)
%   from the provided STA.  Any parameter of the fit can be held fixed or
%   allowed to vary.  A vector of the fit parameters is constructed and
%   passed to the function STA_FIT_ERROR.m  
%
%  AUTHOR: GDF
%  DATE: 2011-06-11

% This information is for debuggin purposes:  sta_fit_error.m is handed a
% long list of unlabled parameters.  The order of these parameter is
% important.  Below is how the parameters should be fed to this function.
% Any change to sta_fit_error, sta_fit_function, or this function should make sure
% that this meaning of the parameter list is not broken.  
%
% PARAMETER FIT INFORMATION 
% 1 center_point_x
% 2 center_point_y
% 3 sd_x
% 4 sd_y
% 5 rotation_angle
% 6 color_weight_a
% 7 color_weight_b
% 8 color weight_c
% 9 x_dim
% 10 y_dim
% 11 num_colors
% 12 fit_surround
% 13 surround_sd_scale
% 14 surround_amp_scale (relative to center_amp)
% 15 scale_one
% 16 scale_two
% 17 tau_one
% 18 tau_two
% 19 n-filters
% 20 frame_number



% ----- PARSE INPUTS -----

p = inputParser;
p.addRequired('sta', @isnumeric);

% optional information to pass (will use if passed)
p.addParamValue('sig_stixels', []);
p.addParamValue('time_course', [], @isnumeric);

% dimensions of STA
p.addParamValue('x_dim', size(sta,2), @isnumeric);
p.addParamValue('y_dim', size(sta,1), @isnumeric);
p.addParamValue('num_colors', size(sta,3), @isnumeric);
p.addParamValue('frame_number', size(sta,4), @isnumeric);

% parameters for obtaining marks
p.addParamValue('mark_params', [], @isstruct);

% parameters for obtaining time course
p.addParamValue('tc_params', [], @isstruct);

p.addParamValue('initial_color_weight_a', [], @isnumeric);
p.addParamValue('initial_color_weight_b', [], @isnumeric);
p.addParamValue('initial_color_weight_c', [], @isnumeric);
p.addParamValue('color_to_normalize', 'peak', @(x)any(1,2,3,'peak'))

% initial values for center
p.addParamValue('initial_center_point_x', [], @isnumeric);
p.addParamValue('initial_center_point_y', [], @isnumeric);
p.addParamValue('initial_center_sd_x', [], @isnumeric);
p.addParamValue('initial_center_sd_y', [], @isnumeric);
p.addParamValue('initial_center_rotation_angle', [], @isnumeric);

% initial values for surround
p.addParamValue('initial_surround_sd_scale', 2, @isnumeric);
p.addParamValue('initial_surround_amp_scale', 0.666/sqrt(2), @isnumeric); % this equation is empirical, assumes surround_amp_scale = 2

% identify what to vary
p.addParamValue('fit_center_point_x', true, @islogical);
p.addParamValue('fit_center_point_y', true, @islogical);
p.addParamValue('fit_center_sd_x', true, @islogical);
p.addParamValue('fit_center_sd_y', true, @islogical);
p.addParamValue('fit_center_rotation_angle', true, @islogical);
p.addParamValue('fit_center_amp_scale', true, @islogical);
p.addParamValue('fit_color_weight_a', true, @islogical);
p.addParamValue('fit_color_weight_b', true, @islogical);
p.addParamValue('fit_color_weight_c', true, @islogical);

p.addParamValue('fit_surround', false, @islogical);
p.addParamValue('fit_surround_sd_scale', false, @islogical);
p.addParamValue('fit_surround_amp_scale', false, @islogical);

% time course information
p.addParamValue('initial_scale_one', [], @isnumeric);
p.addParamValue('initial_scale_two', [], @isnumeric);
p.addParamValue('initial_tau_one', [], @isnumeric);
p.addParamValue('initial_tau_two', [], @isnumeric);
p.addParamValue('initial_n_one_filters', 6, @isnumeric);
p.addParamValue('initial_n_two_filters', 6, @isnumeric);

% time course parameters to vary
p.addParamValue('fit_scale_one', true, @islogical);
p.addParamValue('fit_scale_two', true, @islogical);
p.addParamValue('fit_tau_one', true, @islogical);
p.addParamValue('fit_tau_two', true, @islogical);
p.addParamValue('fit_n_one_filters', false, @islogical);
p.addParamValue('fit_n_two_filters', false, @islogical);

% fit_only_sig stixels
p.addParamValue('fit_sig_stixels_only', false, @islogical);

p.addParamValue('verbose', false, @islogical);

% fiting options
p.addParamValue('optim',{'TolFun',0.001,'Display','off', 'MaxIter', 10000, 'MaxFunEvals', 10000});

% Sig stixels option
p.addParamValue('biggest_blob', true, @islogical);

p.parse(sta, varargin{:});

initial_center_point_x = p.Results.initial_center_point_x;
initial_center_point_y = p.Results.initial_center_point_y;
initial_center_sd_x = p.Results.initial_center_sd_x;
initial_center_sd_y = p.Results.initial_center_sd_y;
initial_center_rotation_angle = p.Results.initial_center_rotation_angle;
initial_color_weight_a = p.Results.initial_color_weight_a;
initial_color_weight_b = p.Results.initial_color_weight_b;
initial_color_weight_c = p.Results.initial_color_weight_c;
initial_surround_sd_scale = p.Results.initial_surround_sd_scale;
initial_surround_amp_scale = p.Results.initial_surround_amp_scale;
initial_scale_one = p.Results.initial_scale_one;
initial_scale_two = p.Results.initial_scale_two;
initial_tau_one = p.Results.initial_tau_one;
initial_tau_two = p.Results.initial_tau_two;
initial_n_one_filters = p.Results.initial_n_one_filters;
initial_n_two_filters = p.Results.initial_n_two_filters;


fit_center_point_x = p.Results.fit_center_point_x;
fit_center_point_y = p.Results.fit_center_point_y;
fit_center_sd_x = p.Results.fit_center_sd_x;
fit_center_sd_y = p.Results.fit_center_sd_y;
fit_center_rotation_angle = p.Results.fit_center_rotation_angle;
fit_color_weight_a = p.Results.fit_color_weight_a;
fit_color_weight_b = p.Results.fit_color_weight_b;
fit_color_weight_c = p.Results.fit_color_weight_c;
fit_surround = p.Results.fit_surround;
fit_surround_sd_scale = p.Results.fit_surround_sd_scale;
fit_surround_amp_scale = p.Results.fit_surround_amp_scale;
fit_scale_one = p.Results.fit_scale_one;
fit_scale_two = p.Results.fit_scale_two;
fit_tau_one = p.Results.fit_tau_one;
fit_tau_two = p.Results.fit_tau_two;
fit_n_one_filters = p.Results.fit_n_one_filters;
fit_n_two_filters = p.Results.fit_n_two_filters;

verbose = p.Results.verbose;
fit_sig_stixels_only = p.Results.fit_sig_stixels_only;

mark_params = p.Results.mark_params;
tc_params = p.Results.tc_params;
x_dim = p.Results.x_dim;
y_dim = p.Results.y_dim;
num_colors = p.Results.num_colors;
frame_number = p.Results.frame_number;
color_to_normalize = p.Results.color_to_normalize;
biggest_blob = p.Results.biggest_blob; 

optim = p.Results.optim;

% check to see if STA is BW
if size(sta, 3) == 1
    fit_color_weight_a = false;
    fit_color_weight_b = false;
    fit_color_weight_c = false;
end


% ---- COMPUTE INITIAL CONDITIONS FOR THE FIT ----



%% GET SIGNIFICANT STIXELS IF NOT PROVIDED

% check to see if sig_stixels (marks) were provided.
% get significant stixels if they were not provided
if isempty(p.Results.sig_stixels)
    if isempty(mark_params)
        sig_stixels = significant_stixels(sta);
    else
        sig_stixels =significant_stixels(sta, mark_params);
    end
else
    sig_stixels = p.Results.sig_stixels;
end

% make sure there are significant stixels, otherwise report error
if sum(sum(sig_stixels)) == 0
    warning('no significant stixels were found or provided for STA, lower threshold or provide sig stixels')
    fit_info = [];
    return
end


% get matrix subscripts to these pixels for the eigenvalue calculation
if biggest_blob
    biggestBlob_stixs = ExtractNLargestBlobs(full(sig_stixels), 1);
    [matrix_subscript_i_eig, matrix_subscript_j_eig] = find(biggestBlob_stixs);
    matrix_subscripts_eig = [matrix_subscript_i_eig, matrix_subscript_j_eig];
else
    [matrix_subscript_i_eig, matrix_subscript_j_eig] = find(sig_stixels);
    matrix_subscripts_eig = [matrix_subscript_i_eig, matrix_subscript_j_eig];
end


% get matrix subscripts of all sig stixels

% [matrix_subscript_i, matrix_subscript_j] = find(sig_stixels);
% matrix_subscripts = [matrix_subscript_i, matrix_subscript_j];




% GET TIME COURSE INFORMATION

% compute the time course from the sig_stixels (marks), unless time course
% was provided
if isempty(p.Results.time_course)
    if isempty(tc_params)
       time_course = time_course_from_sta(sta, sig_stixels);
    else
       time_course = time_course_from_sta(sta, sig_stixels, tc_params);
    end
else
    time_course = p.Results.time_course;
end

%% find the peak and trough of TC

% collapse TC to one color if RGB
if size(time_course, 2) > 1;
    summed_tc = sum(time_course,2); % collapse RGB tc
else
    summed_tc = time_course; % don't collapse a BW tc
end

% matlab indexes the time course such that the highest index is the frame
% closest to the spike time. The fit function assumes the opposite. 
% To compensate the min(max) frame is subtracted from the total frame
% number
[max_tc_val, max_frame] = max(summed_tc);
[min_tc_val, min_frame] = min(summed_tc);
max_frame = frame_number - max_frame;
min_frame = frame_number - min_frame;

% Sort parameters appropriately for ON and OFF cells
if min_frame < max_frame % true for off cells
    peak_frame = min_frame;
    trough_frame = max_frame;
    peak_scale = min_tc_val;
    trough_scale = max_tc_val;
else                     % true for on cells
    peak_frame = max_frame;
    trough_frame = min_frame;
    peak_scale = max_tc_val;
    trough_scale = min_tc_val;
end






%% SET INITIAL TEMPORAL FIT VALUES
if isempty(initial_scale_one)
    initial_scale_one = peak_scale * 1.05; % this value is empirical (from monkey)
end
    
if isempty(initial_scale_two)
    initial_scale_two = trough_scale * 0.9; % this value is empiricle (from monkey)

end

if isempty(initial_tau_one)
%     initial_tau_one = peak_frame *1.05;
        initial_tau_one = peak_frame *10;

end

% initial tau_two depends a bit on  trough_scale
if isempty(initial_tau_two)
    if abs(trough_scale/peak_scale) > 0.15
%         initial_tau_two = trough_frame;
                initial_tau_two = trough_frame*10;

    else
%         initial_tau_two = trough_frame * 1.25;
                initial_tau_two = trough_frame * 10;

    end
end


%% GET INITIAL SPATIAL PARAMETERS

% if surround isn't to be fit, then set these parameter to zero
if ~fit_surround
   initial_surround_sd_scale = 0;
   initial_surround_amp_scale = 0;
end

% Get initial center points
if any(isempty([initial_center_point_x, initial_center_point_y]))
    % compute the initial center of spatial fit
    if size(matrix_subscripts_eig,1) > 1
        tmp = centroid(matrix_subscripts_eig);
        initial_center_point_x = tmp(2);
        initial_center_point_y = tmp(1);
    else
        initial_center_point_x = matrix_subscripts_eig(2);
        initial_center_point_y = matrix_subscripts_eig(1);
    end
end
  
% Get initial sigmas (x and y) and rotation angle
if size(matrix_subscripts_eig,1) >= 3
    % base initial conds on PCA permformed on marks
    [correlation_matrix, ~, eigvalues] = princomp(matrix_subscripts_eig);

    if isempty(initial_center_rotation_angle)
        initial_center_rotation_angle = pi/2 * correlation_matrix(2);
    end

    if isempty(initial_center_sd_x)
        initial_center_sd_x = sqrt(eigvalues(1)) * 1.2;
        if initial_center_sd_x < 0.5
            initial_center_sd_x = 0.5;
        end
    end
    
    if isempty(initial_center_sd_y)
        initial_center_sd_y = sqrt(eigvalues(2)) * 1.2;
        if initial_center_sd_y < 0.5
            initial_center_sd_y = 0.5;
        end
    end
else
    % if there aren't many marks, start with small circular RF.
    if isempty(initial_center_rotation_angle)
        initial_center_rotation_angle = 0;
    end
    
    if isempty(initial_center_sd_x)
        initial_center_sd_x = 1;
    end
    
    if isempty(initial_center_sd_y)
        initial_center_sd_y = 1;
    end
    
end


%% INITIAL COLORS

% identify the dominant color
if size(sta, 3) == 3;
    [~, peak_color] = max(abs(time_course(frame_number - peak_frame,:)));


    switch color_to_normalize
        case 'peak' 

            % get initial peak color triplet and normalize to unit length
            peak_color_triplet = abs(time_course((frame_number - peak_frame),:));
            peak_color_triplet = peak_color_triplet ./ norm(peak_color_triplet);
            
            initial_color_weight_a = peak_color_triplet(1);
            initial_color_weight_b = peak_color_triplet(2);
            initial_color_weight_c = peak_color_triplet(3);
            
    end
    
elseif size(sta, 3) == 1; % for Black&White STAs

    initial_color_weight_a = 1;
    initial_color_weight_b = 1;
    initial_color_weight_c = 1;
else
    error('dimensions of STA color not recognized')
end




%% --- Sort out which parameters to vary and which to hold fixed ---

% a list of all the input parameters to be passed to fit function

input_params = [...
    initial_center_point_x,...
    initial_center_point_y,...
    initial_center_sd_x,...
    initial_center_sd_y,...
    initial_center_rotation_angle,...
    initial_color_weight_a,...
    initial_color_weight_b,...
    initial_color_weight_c,...
    x_dim,...
    y_dim,....
    num_colors,...
    fit_surround,...
    initial_surround_sd_scale,...
    initial_surround_amp_scale,...
    initial_scale_one,...
    initial_scale_two,...
    initial_tau_one,...
    initial_tau_two,...
    initial_n_one_filters,...
    frame_number,...
    initial_n_two_filters,...
    ];
    
% make a logical matrix that identifies which parameters are free
fit_list = [...
    fit_center_point_x,...
    fit_center_point_y,...
    fit_center_sd_x,...
    fit_center_sd_y,...
    fit_center_rotation_angle,...
    fit_color_weight_a,...
    fit_color_weight_b,...
    fit_color_weight_c,...
    false,... x_dim is not fit
    false,... y_dim is not fit
    false,... num_colors is not fit
    false,... even when fit_surround is true, we do not want its value changed
    fit_surround_sd_scale,...
    fit_surround_amp_scale,...
    fit_scale_one,...
    fit_scale_two,...
    fit_tau_one,...
    fit_tau_two,...
    fit_n_one_filters,...
    false,... % frame number is not fit
    fit_n_two_filters,...
    ];   


% logical indices that will identify which parameters to vary vs. fix.
fit_indices = find(fit_list);
fixed_indices = find(~fit_list);




%% --- FIT THE STA ---

input_params = double(input_params);

[final_fit_params, fval] = fminsearch(@(fit_params)sta_fit_error(sta,...
                              fit_params, input_params(fixed_indices),...
                              fit_indices, fixed_indices,...
                              verbose, mark_params),...
                              input_params(fit_indices),...
                              optimset(optim{:}));
                          
final_params = input_params;
final_params(fit_indices) = final_fit_params;
sta_fit = sta_fit_function(final_params);

                          
%% parse the output and organize into a structure

temp_pointer = 1;

% center_x
if fit_list(1)
    fit_info.center_point_x = final_fit_params(temp_pointer);
    fit_info.fit_center_point_x = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.center_point_x = initial_center_point_x;
    fit_info.fit_center_point_x = false;
end

% center_y;
if fit_list(2)
    fit_info.center_point_y = final_fit_params(temp_pointer);
    fit_info.fit_center_point_y = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.center_point_y = initial_center_point_y;
    fit_info.fit_center_point_y = false;
end

% center_sd_x
if fit_list(3)
    fit_info.center_sd_x = final_fit_params(temp_pointer);
    fit_info.fit_center_sd_x = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.center_sd_x = initial_center_sd_x;
    fit_info.fit_center_sd_x = false;
end

% center_sd_y
if fit_list(4)
    fit_info.center_sd_y = final_fit_params(temp_pointer);
    fit_info.fit_center_sd_y = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.center_sd_y = initial_center_sd_y;
    fit_info.fit_center_sd_y = false;
end

% center_rotation_angle
if fit_list(5)
    fit_info.center_rotation_angle = final_fit_params(temp_pointer);
    fit_info.fit_center_rotation_angle = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.center_rotation_angle = initial_center_rotation_angle;
    fit_info.fit_center_rotation_angle = false;
end

% color_a
if fit_list(6)
    fit_info.color_weight_a = final_fit_params(temp_pointer);
    fit_info.fit_color_weight_a = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.color_weight_a = initial_color_weight_a;
    fit_info.fit_color_weight_a = false;
end

% color_b
if fit_list(7)
    fit_info.color_weight_b = final_fit_params(temp_pointer);
    fit_info.fit_color_weight_b = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.color_weight_b = initial_color_weight_b;
    fit_info.fit_color_weight_b = false;
end

% color_c
if fit_list(8)
    fit_info.color_weight_c = final_fit_params(temp_pointer);
    fit_info.fit_color_weight_c = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.color_weight_c = initial_color_weight_c;
    fit_info.fit_color_weight_c = false;
end

fit_info.y_dim = y_dim;
fit_info.x_dim = x_dim;
fit_info.num_colors = num_colors;
fit_info.color_to_normalize = color_to_normalize;

if fit_list(13)
    fit_info.surround_sd_scale = final_fit_params(temp_pointer);
    fit_info.fit_surround_sd_scale = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.surround_sd_scale = initial_surround_sd_scale;
    fit_info.fit_surround_sd_scale = false;
end

if fit_list(14)
    fit_info.surround_amp_scale = final_fit_params(temp_pointer);
    fit_info.fit_surround_amp_scale = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.surround_amp_scale = initial_surround_amp_scale;
    fit_info.fit_surround_amp_scale = false;
end

% tc scale one
if fit_list(15)
    fit_info.scale_one = final_fit_params(temp_pointer);
    fit_info.fit_scale_one = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.scale_one = initial_scale_one;
    fit_info.fit_scale_one = false;
end

% tc scale two
if fit_list(16)
    fit_info.scale_two = final_fit_params(temp_pointer);
    fit_info.fit_scale_two = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.scale_two = initial_scale_two;
    fit_info.fit_scale_two = false;
end

% tc tau one
if fit_list(17)
    fit_info.tau_one = final_fit_params(temp_pointer);
    fit_info.fit_tau_one = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.tau_one = initial_tau_one;
    fit_info.fit_tau_one = false;
end

% tc tau two
if fit_list(18)
    fit_info.tau_two = final_fit_params(temp_pointer);
    fit_info.fit_tau_two = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.tau_two = initial_tau_two;
    fit_info.fit_tau_two = false;
end

% n one filters
if fit_list(19)
    fit_info.n_one_filters = final_fit_params(temp_pointer);
    fit_info.fit_n_one_filters = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.n_one_filters = initial_n_one_filters;
    fit_info.fit_n_one_filters = false;
end

% frame number
if fit_list(20)
    fit_info.frame_number = final_fit_params(temp_pointer);
    fit_info.frame_number = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.frame_number= frame_number;
    fit_info.fit_frame_number = false;
end

% n two filters
if fit_list(21)
    fit_info.n_two_filters = final_fit_params(temp_pointer);
    fit_info.fit_n_two_filters = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.n_two_filters = initial_n_two_filters;
    fit_info.fit_n_two_filters = false;
end


fit_info.fit_surround = fit_surround;
fit_info.frame_number = frame_number;
fit_info.rmse = fval;
fit_info.initial_params = input_params;
    fit_info.fit_indices = fit_indices;
    fit_info.fixed_indices = fixed_indices;
    fit_info.fit_params = final_fit_params;
    fit_info.fixed_params = input_params(fixed_indices);

