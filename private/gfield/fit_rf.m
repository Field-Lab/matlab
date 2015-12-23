function [fit_info, initial_params] = fit_sta(sta, varargin)
% fit_info = fit_sta(sta, varargin)
%
% hand the function an STA and it will fit the movie with a differenc of
% Gaussian function in space and the product of two filters in time
%
% 
% A list of parameters, some free to vary, and some fixed are passed
% to fit functions.  To reduce compute time, these fit functions are
% dumb and must be passed the correct numbers in the correct order.
% Therefore, a parameter list is generated from the user provided inputs
% that contains all the fit parameters and this must be correctly indexed
% when passing the information in this parameter list to the lower level
% functions.  If it is not correctly index, you might think you are fitting
% the position of the RF, when in fact you are fitting its width.  Thus,
% for debugging purposes, a master list of the all the parameter and their
% index in the parameter list (called all_params) is provided here.  Note
% every value passed to the fit is either an scalar or a logical (this
% allows the low level functions to be as dumb, and thus as fast, as
% possible).
%
% FIT INFORMATION 
% 1 center_point_x
% 2 center_point_y
% 3 sd_x
% 4 sd_y
% 5 center_amp
% 6 rotation_angle
% 7 color_weight_a
% 8 color_weight_b
% 9 color weight_c
% 10 x_dim
% 11 y_dim
% 12 num_colors
% 13 color_to_normalize
% 14 surround_sd_scale
% 15 surround_amp_scale (relative to center_amp)


p = inputParser;
p.addRequired('sta', @isnumeric);

% dimensions of STA
p.addParamValue('x_dim', size(sta,2), @isnumeric);
p.addParamValue('y_dim', size(sta,1), @isnumeric);
p.addParamValue('num_colors', size(sta,3), @isnumeric);

% parameters for obtaining marks
p.addParamValue('mark_params', [], @isstruct);

p.addParamValue('initial_color_weight_a', 1, @isnumeric);
p.addParamValue('initial_color_weight_b', 1, @isnumeric);
p.addParamValue('initial_color_weight_c', 1, @isnumeric);
p.addParamValue('color_to_normalize', 'peak', @(x)any(1,2,3,'peak'))

% initial values for center
p.addParamValue('initial_center_point_x', [], @isnumeric);
p.addParamValue('initial_center_point_y', [], @isnumeric);
p.addParamValue('initial_center_sd_x', [], @isnumeric);
p.addParamValue('initial_center_sd_y', [], @isnumeric);
p.addParamValue('initial_center_rotation_angle', [], @isnumeric);
p.addParamValue('initial_center_amp_scale', [], @isnumeric)

% initial values for surround
p.addParamValue('initial_surround_sd_scale', 2, @isnumeric);
p.addParamValue('initial_surround_amp_scale', 0.2, @isnumeric);

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

p.parse(sta, varargin{:});

initial_center_point_x = p.Results.initial_center_point_x;
initial_center_point_y = p.Results.initial_center_point_y;
initial_center_sd_x = p.Results.initial_center_sd_x;
initial_center_sd_y = p.Results.initial_center_sd_y;
initial_center_amp_scale = p.Results.initial_center_amp_scale;
initial_center_rotation_angle = p.Results.initial_center_rotation_angle;
initial_color_weight_a = p.Results.initial_color_weight_a;
initial_color_weight_b = p.Results.initial_color_weight_b;
initial_color_weight_c = p.Results.initial_color_weight_c;
initial_surround_sd_scale = p.Results.initial_surround_sd_scale;
initial_surround_amp_scale = p.Results.initial_surround_amp_scale;

fit_center_point_x = p.Results.fit_center_point_x;
fit_center_point_y = p.Results.fit_center_point_y;
fit_center_sd_x = p.Results.fit_center_sd_x;
fit_center_sd_y = p.Results.fit_center_sd_y;
fit_center_rotation_angle = p.Results.fit_center_rotation_angle;
fit_center_amp_scale = p.Results.fit_center_amp_scale;
fit_color_weight_a = p.Results.fit_color_weight_a;
fit_color_weight_b = p.Results.fit_color_weight_b;
fit_color_weight_c = p.Results.fit_color_weight_c;
fit_surround = p.Results.fit_surround;
fit_surround_sd_scale = p.Results.fit_surround_sd_scale;
fit_surround_amp_scale = p.Results.fit_surround_amp_scale;

mark_params = p.Results.mark_params;
x_dim = p.Results.x_dim;
y_dim = p.Results.y_dim;
num_colors = p.Results.num_colors;
color_to_normalize = p.Results.color_to_normalize;

robust_percent = 0.75;


%% GET SIGNIFICANT STIXELS

% get significant stixels to initialize parameters
if isempty(mark_params)
    sig_stixels = significant_stixels(sta);
else
    sig_stixels = significant_stixels(sta, mark_params);
end
% get matrix subscripts to these pixels
[matrix_subscript_i, matrix_subscript_j] = find(sig_stixels);
matrix_subscripts = [matrix_subscript_i, matrix_subscript_j];


% find the peak frame
time_course = time_course_from_sta(sta, sig_stixels);
if size(time_course, 2) > 1;
    summed_time_course = sum(time_course,2);
    [~, peak_frame] = max(abs(summed_time_course));
end


%% GET INITIAL SPATIAL PARAMETERS

% center_point
if any(isempty([initial_center_point_x, initial_center_point_y]))
    % compute the initial center of spatial fit
    if size(matrix_subscripts,1) > 1
        tmp = centroid(matrix_subscripts);
        initial_center_point_x = tmp(2);
        initial_center_point_y = tmp(1);
    else
        initial_center_point_x = matrix_subscripts(2);
        initial_center_point_y = matrix_subscripts(1);
    end
end
  


% center_sd_scale
if any(isempty([initial_center_sd_x, initial_center_sd_y]))
    if size(matrix_subscripts,1) > 1;
        distance_list = ipdm(matrix_subscripts);
        distance_list = triu(distance_list);
        [~, ~, distance_list] = find(distance_list);
        initial_center_sd_scale = robust_mean(distance_list, robust_percent);
        initial_center_sd_x = initial_center_sd_scale;
        initial_center_sd_y = initial_center_sd_scale;
    else
        initial_center_sd_x = 1;
        initial_center_sd_y = 1;
    end
end



% center_rotation_angle
if isempty(initial_center_rotation_angle)
    initial_center_rotation_angle = 0.0;
end



% center_amp_scale
if isempty(initial_center_amp_scale)
    initial_center_amp_scale = ext(reshape(sta, 1,[]));
end

% color weights
switch color_to_normalize
    case 'peak' 
        tmp = sta(round(initial_center_point_y), round(initial_center_point_x), :, peak_frame);
        [~, extrema_index] = max(abs(squeeze(tmp)));
        
        switch extrema_index
            case 1
                initial_color_weight_a = 1;
                initial_color_weight_b = initial_color_weight_b ./ initial_color_weight_a;
                initial_color_weight_c = initial_color_weight_c ./ initial_color_weight_a;
                fit_color_weight_a = false;
            case 2
                initial_color_weight_b = 1;
                initial_color_weight_a = initial_color_weight_a ./ initial_color_weight_b;
                initial_color_weight_c = initial_color_weight_c ./ initial_color_weight_b;
                fit_color_weight_b = false;
            case 3 
                initial_color_weight_c = 1;
                initial_color_weight_a = initial_color_weight_a ./ initial_color_weight_c;
                initial_color_weight_b = initial_color_weight_b ./ initial_color_weight_c;
                fit_color_weight_c = false;
        end
        
    case 1
        initial_color_weight_a = 1;
        initial_color_weight_b = initial_color_weight_b ./ initial_color_weight_a;
        initial_color_weight_c = initial_color_weight_c ./ initial_color_weight_a;
        fit_color_weight_a = false;
    case 2
        initial_color_weight_b = 1;
        initial_color_weight_a = initial_color_weight_a ./ initial_color_weight_b;
        initial_color_weight_c = initial_color_weight_c ./ initial_color_weight_b;
        fit_color_weight_b = false;
    case 3 
        initial_color_weight_c = 1;
        initial_color_weight_a = initial_color_weight_a ./ initial_color_weight_c;
        initial_color_weight_b = initial_color_weight_b ./ initial_color_weight_c;
        fit_color_weight_c = false;
end


%%


%% Sort out which parameters to vary and which to hold fixed


input_params = [...
    initial_center_point_x,...
    initial_center_point_y,...
    initial_center_sd_x,...
    initial_center_sd_y,...
    initial_center_amp_scale,...
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
    ];
    
% make a logical matrix that identifies which parameters are free
fit_list = [...
    fit_center_point_x,...
    fit_center_point_y,...
    fit_center_sd_x,...
    fit_center_sd_y,...
    fit_center_amp_scale,...
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
    ];   


fit_indices = find(fit_list);
fixed_indices = find(~fit_list);


%% FIT THE STA

input_params = double(input_params);

final_fit_params = fminsearch(@(fit_params)sta_fit_error(sta,...
                              fit_params, input_params(fixed_indices),...
                              fit_indices, fixed_indices),...
                              input_params(fit_indices)); 
                          

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

% center_amp_scale
if fit_list(5)
    fit_info.center_amp_scale = final_fit_params(temp_pointer);
    fit_info.fit_center_amp_scale = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.center_amp_scale = initial_center_amp_scale;
    fit_info.fit_center_amp_scale = false;
end

% center_rotation_angle
if fit_list(6)
    fit_info.center_rotation_angle = final_fit_params(temp_pointer);
    fit_info.fit_center_rotation_angle = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.ccenter_rotation_angle = initial_center_rotation_angle;
    fit_info.center_rotation_angle = false;
end

% color_a
if fit_list(7)
    fit_info.color_weight_a = final_fit_params(temp_pointer);
    fit_info.fit_color_weight_a = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.color_weight_a = initial_color_weight_a;
    fit_info.fit_color_weight_a = false;
end

% color_b
if fit_list(8)
    fit_info.color_weight_b = final_fit_params(temp_pointer);
    fit_info.fit_color_weight_b = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.color_weight_b = initial_color_weight_b;
    fit_info.fit_color_weight_b = false;
end

% color_c
if fit_list(9)
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
fit_info.fit_surround = fit_surround;

if fit_list(14)
    fit_info.surround_sd_scale = final_fit_params(temp_pointer);
    fit_info.fit_surround_sd_scale = true;
    temp_pointer = temp_pointer + 1;
else
    fit_info.surround_sd_scale = initial_surround_sd_scale;
    fit_info.fit_surround_sd_scale = false;
end

if fit_list(15)
    fit_info.surround_amp_scale = final_fit_params(temp_pointer);
    fit_info.fit_surround_amp_scale = true;
    temp_pointer = temp_pointer +1;
else
    fit_info.surround_amp_scale = initial_surround_amp_scale;
    fit_info.fit_surround_amp_scale = fit_surround_amp_scale;
end

initial_params = input_params;



    
    








   

