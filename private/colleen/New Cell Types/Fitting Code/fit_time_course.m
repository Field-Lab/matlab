function [tc_fit, final_params] = fit_time_course(time_course, varargin)
%
% The function fits a BW time course vector.  time_course_fit_function and 
% tc_fit_error are functions that are called.  The fit is a difference of
% nth-order low pass filters .
%
% Usage: [tc_fit, final_params] = fit_time_course(time_course, varargin)
%
% inputs: 
%    time_course               BW time course column vector
% 
% optional inputs
%   initial parameters for the fit can be passed as well as instructions on
%   which parameters to vary and which to keep fix. 
%   See below in the inputParse definitions for more information about
%   defaults
%   
%   optim structure can be passed to control quality of the fit.
%
% outputs: 
%   tc_fit                  time course fit
%   final_params            the final fit parameters
%
%
% Created: GDF 5/2014
%

p = inputParser;
p.addRequired('time_course', @isnumeric);

% time course information
p.addParamValue('initial_scale_one', [], @isnumeric);
p.addParamValue('initial_scale_two', [], @isnumeric);
p.addParamValue('initial_tau_one', [], @isnumeric);
p.addParamValue('initial_tau_two', [], @isnumeric);
p.addParamValue('initial_n_filters', 10, @isnumeric);

% time course parameters to vary
p.addParamValue('fit_scale_one', true, @islogical);
p.addParamValue('fit_scale_two', true, @islogical);
p.addParamValue('fit_tau_one', true, @islogical);
p.addParamValue('fit_tau_two', true, @islogical);
p.addParamValue('fit_n_filters', true, @islogical);

p.addParamValue('verbose', false, @islogical);

% fiting options
p.addParamValue('optim',{'TolFun',0.0001,'Display','off', 'MaxIter', 100000, 'MaxFunEvals', 100000});

p.parse(time_course, varargin{:});

frame_number = length(time_course);


initial_scale_one = p.Results.initial_scale_one;
initial_scale_two = p.Results.initial_scale_two;
initial_tau_one = p.Results.initial_tau_one;
initial_tau_two = p.Results.initial_tau_two;
initial_n_filters = p.Results.initial_n_filters;

fit_scale_one = p.Results.fit_scale_one;
fit_scale_two = p.Results.fit_scale_two;
fit_tau_one = p.Results.fit_tau_one;
fit_tau_two = p.Results.fit_tau_two;
fit_n_filters = p.Results.fit_n_filters;

verbose = p.Results.verbose;

optim = p.Results.optim;



%% find the peak and trough of TC

% matlab indexes the time course such that the highest index is the frame
% closest to the spike time. The fit function assumes the opposite. 
% To compensate the min(max) frame is subtracted from the total frame
% number
[max_tc_val, max_frame] = max(time_course);
[min_tc_val, min_frame] = min(time_course);
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
    initial_scale_two = trough_scale * 0.95; % this value is empiricle (from monkey)
end

if isempty(initial_tau_one)
    initial_tau_one = peak_frame * 1.1;
end

% initial tau_two depends a bit on trough_scale
if isempty(initial_tau_two)
    if abs(trough_scale/peak_scale) > 0.15
        initial_tau_two = trough_frame;
    else
        initial_tau_two = trough_frame * 1.25;
    end
end

%% --- Sort out which parameters to vary and which to hold fixed ---

% a list of all the input parameters to be passed to fit function

input_params = [...
    initial_scale_one,...
    initial_scale_two,...
    initial_tau_one,...
    initial_tau_two,...
    initial_n_filters,...
    frame_number,...
    ];
    
% make a logical matrix that identifies which parameters are free
fit_list = [...
    fit_scale_one,...
    fit_scale_two,...
    fit_tau_one,...
    fit_tau_two,...
    fit_n_filters,...
    false,... % frame number is not fit
    ];   


% logical indices that will identify which parameters to vary vs. fix.
fit_indices = find(fit_list);
fixed_indices = find(~fit_list);



%% --- FIT THE time course ---

input_params = double(input_params);

[final_fit_params, fval] = fminsearch(@(fit_params)tc_fit_error(time_course,...
                              fit_params, input_params(fixed_indices),...
                              fit_indices, fixed_indices,...
                              verbose),...
                              input_params(fit_indices),...
                              optimset(optim{:}));
                          
final_params = input_params;
final_params(fit_indices) = final_fit_params;
tc_fit = time_course_fit_function(final_params);

