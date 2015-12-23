function fit_info = fit_sta_sequence(sta, varargin)
% FIT_STA_SEQUENCE     compute the fit for sta by sequentially 
%                      fitting temporal, center and surround.
%
% usage:  fit_info = fit_sta_sequence(sta, varargin)
%
% INPUTS:
%   sta              a spatial-temporal-chromatic STA
%
% OPTIONAL INPUTS:
%   fit_temporal       true
%   fit_center         true
%   fit_surround       true
%   fit_instructions   []
%   verbose            false      
%   
%
%
% OUPUTS
%   fit_info        A structure that contains all the information about
%                       the fit, including initial conditions
%
%
%
% 2013-04 xyao
%



p = inputParser;

p.addRequired('sta', @isnumeric);

p.addParamValue('fit_temporal', true, @islogical);
p.addParamValue('fit_center', true, @islogical);
p.addParamValue('fit_surround', true, @islogical);
p.addParamValue('fit_instructions', [], @isstruct);

p.addParamValue('verbose', false, @islogical);


p.parse(sta, varargin{:});

fit_temporal = p.Results.fit_temporal;
fit_center = p.Results.fit_center;
fit_surround = p.Results.fit_surround;
verbose = p.Results.verbose;
fit_instructions = p.Results.fit_instructions;



%% BODY OF THE FUNCTION


% set thresh for pick out sig stixel.
mark_params.thresh = 4.0;


% fit sta sequentially 
temp_sta = sta;
    
    % hold spatial params fixed
    if fit_temporal
        clear params
        params.fit_center_point_x = false;
        params.fit_center_point_y = false;
        params.fit_center_sd_x = false;
        params.fit_center_sd_y = false;
        params.fit_center_rotation_angle = false;
        params.fit_center_amp_scale = false;
        params.fit_surround = false;
        params.fit_surround_sd_scale = false;
        params.fit_surround_amp_scale = false;
        params.verbose = verbose;

        params.mark_params = mark_params;

        % fit an STA with matlab code.
        if isempty(fit_instructions)
            temp_fit = fit_sta(temp_sta, params);
        else
            temp_fit = fit_sta(temp_sta, params, fit_instructions);
        end
        
    end
    
    % use previous fit for temporal params and hold them fixed -- fit center
    if fit_center
        clear params
        if fit_temporal
            % set params
            params.initial_scale_one = temp_fit.scale_one;
            params.initial_scale_two = temp_fit.scale_two;
            params.initial_tau_one = temp_fit.tau_one;
            params.initial_tau_two = temp_fit.tau_two;
            params.initial_n_filters = temp_fit.n_filters;
        end
        
            % fix params
            params.fit_scale_one = false;
            params.fit_scale_two = false;
            params.fit_tau_one = false;
            params.fit_tau_two = false;
            params.fit_n_filters = false;
            params.verbose = verbose;

            % set spatial surround to 0 and don't fit
            params.initial_surround_sd_scale = 0;
            params.initial_surround_amp_scale = 0;
            params.fit_surround = false;
            params.fit_surround_sd_scale = false;
            params.fit_surround_amp_scale = false;

            params.mark_params = mark_params;
    
            if isempty(fit_instructions)
                temp_fit = fit_sta(temp_sta, params);
            else
                temp_fit = fit_sta(temp_sta, params, fit_instructions);
            end
    end
    
    % Introduce surrround to fit while keeping everything else constant 
    if fit_surround
        clear params
        if  fit_center
            % set center parameters 
            params.initial_center_point_x = temp_fit.center_point_x;
            params.initial_center_point_y = temp_fit.center_point_y;
            params.initial_center_sd_x = temp_fit.center_sd_x;
            params.initial_center_sd_y = temp_fit.center_sd_y;
            params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
            
        end
        
        if fit_temporal
            %set temporal parameters 
            params.initial_scale_one = temp_fit.scale_one;
            params.initial_scale_two = temp_fit.scale_two;
            params.initial_tau_one = temp_fit.tau_one;
            params.initial_tau_two = temp_fit.tau_two;
            params.initial_n_filters = temp_fit.n_filters;
        end
        
        % hold temporal and center params fixed
        params.fit_center_point_x = false;
        params.fit_center_point_y = false;
        params.fit_center_sd_x = false;
        params.fit_center_sd_y = false;
        params.fit_center_rotation_angle = false;
        params.fit_center_amp_scale = false;
            
        params.fit_scale_one = false;
        params.fit_scale_two = false;
        params.fit_tau_one = false;
        params.fit_tau_two = false;
        params.fit_n_filters = false;
            
        % fit surround
        params.initial_surround_sd_scale = 1.5;
        params.initial_surround_amp_scale = 0.1;
        params.fit_surround = true;
        params.fit_surround_sd_scale = true;
        params.fit_surround_amp_scale = true;
        params.verbose = verbose;

        params.mark_params = mark_params;

            if isempty(fit_instructions)
                temp_fit = fit_sta(temp_sta, params);
            else
                temp_fit = fit_sta(temp_sta, params, fit_instructions);
            end
    end
 
  % Use these above fit results as initial conditions, then refit temporal
    if fit_temporal  
        clear params
        % set fit params
        params.fit_center_point_x = false;
        params.fit_center_point_y = false;
        params.fit_center_sd_x = false;
        params.fit_center_sd_y = false;
        params.fit_center_rotation_angle = false;
        params.fit_center_amp_scale = false;
        params.fit_surround = true;
        params.fit_surround_sd_scale = false;
        params.fit_surround_amp_scale = false;
        
        % set initial condition
        params.initial_scale_one = temp_fit.scale_one;
        params.initial_scale_two = temp_fit.scale_two;
        params.initial_tau_one = temp_fit.tau_one;
        params.initial_tau_two = temp_fit.tau_two;
        params.initial_n_filters = temp_fit.n_filters;
        
        if fit_center
            params.initial_center_point_x = temp_fit.center_point_x;
            params.initial_center_point_y = temp_fit.center_point_y;
            params.initial_center_sd_x = temp_fit.center_sd_x;
            params.initial_center_sd_y = temp_fit.center_sd_y;
            params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
        end
        
        if fit_surround
            params.initial_surround_sd_scale = temp_fit.surround_sd_scale;
            params.initial_surround_amp_scale = temp_fit.surround_amp_scale;
        end
        
        params.verbose = verbose;
        params.mark_params = mark_params;
        
        % fit an STA with matlab code.
            if isempty(fit_instructions)
                temp_fit = fit_sta(temp_sta, params);
            else
                temp_fit = fit_sta(temp_sta, params, fit_instructions);
            end
    end
    
  % Use these above fit results as initial conditions and allow the entire fit to vary
    clear params
    params.fit_center_point_x = false;
    params.fit_center_point_y = false;
    params.fit_center_sd_x = false;
    params.fit_center_sd_y = false;
    params.fit_center_rotation_angle = false;
        
    params.fit_scale_one = false;
    params.fit_scale_two = false;
    params.fit_tau_one = false;
    params.fit_tau_two = false;
    params.fit_n_filters = false;
    
    if fit_center
        % set center parameters
        params.initial_center_point_x = temp_fit.center_point_x;
        params.initial_center_point_y = temp_fit.center_point_y;
        params.initial_center_sd_x = temp_fit.center_sd_x;
        params.initial_center_sd_y = temp_fit.center_sd_y;
        params.initial_center_rotation_angle = temp_fit.center_rotation_angle;
        params.fit_center_point_x = true;
        params.fit_center_point_y = true;
        params.fit_center_sd_x = true;
        params.fit_center_sd_y = true;
        params.fit_center_rotation_angle = true;
        params.fit_center_amp_scale = true;
    end
    
    if fit_temporal
        % set temporal parameters 
        params.initial_scale_one = temp_fit.scale_one;
        params.initial_scale_two = temp_fit.scale_two;
        params.initial_tau_one = temp_fit.tau_one;
        params.initial_tau_two = temp_fit.tau_two;
        params.initial_n_filters = temp_fit.n_filters;
        params.fit_scale_one = true;
        params.fit_scale_two = true;
        params.fit_tau_one = true;
        params.fit_tau_two = true;
        params.fit_n_filters = true;
    end
    
    if fit_surround
        % set surround parameters
       params.initial_surround_sd_scale = temp_fit.surround_sd_scale;
       params.initial_surround_amp_scale = temp_fit.surround_amp_scale;
       params.fit_surround = true;
       params.fit_surround_sd_scale = true;
       params.fit_surround_amp_scale = true;
    end
    
       params.verbose = verbose;

       params.mark_params = mark_params;
      
            if isempty(fit_instructions)
                temp_fit = fit_sta(temp_sta, params);
            else
                temp_fit = fit_sta(temp_sta, params, fit_instructions);
            end
  
       fit_info = temp_fit;
