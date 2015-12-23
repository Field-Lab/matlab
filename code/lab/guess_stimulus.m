function stimulus_out = guess_stimulus(stimulus, varargin)
% guess_stimulus     Guess what the stimulus parameters were
%
% usage:  stimulus_out = guess_stimulus(stimulus, varargin)
%
% arguments: stimulus - stimulus struct
%            varargin - struct or list of optional parameters (see below)
%
% outputs:   stimulus_out - stimulus struct with fields added
%
%
% optional params, their default values, and what they specify:
%
% monitor_x       	640         monitor width in pixels
% monitor_y       	320         monitor height in pixels
% monitor_refresh   120         monitor refresh rate, in Hz
%
%
% 2009-04 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('monitor_x',640);
p.addParamValue('monitor_y',480);
p.addParamValue('monitor_refresh',120);
p.addParamValue('display_type', 'crt')

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% initialize with input stimulus
stimulus_out = stimulus;

switch params.display_type
    
    case 'crt'
        % choose defaults based on monitor size
        if all([params.monitor_x params.monitor_y] == [640 480])

            % initialize default struct
            default = struct;

            % choose default values based on aspect ratio
            if stimulus.field_height == stimulus.field_width
                % square
                default.x_start = 160;
                default.x_end = 480;
                default.y_start = 80;
                default.y_end = 400;
            else
                % rectangle
                default.x_start = 0;
                default.x_end = 640;
                default.y_start = 80;
                default.y_end = 400;
            end

        else
            error('I don''t know the default stimulus parameters for a monitor of size %d x %d pixels.',params.monitor_x,params.monitor_y)
        end

        % enter missing parameters using the defaults, being careful not to replace any existing fields

        if ~isfield(stimulus_out,'stixel_height')
            stimulus_out.stixel_height = (default.y_end - default.y_start)/stimulus.field_height;
        end

        if ~isfield(stimulus_out,'stixel_width')
            stimulus_out.stixel_width = (default.x_end - default.x_start)/stimulus.field_width;
        end

        if ~isfield(stimulus_out,'x_start'),stimulus_out.x_start = default.x_start; end

        if ~isfield(stimulus_out,'x_end'), stimulus_out.x_end = default.x_end; end

        if ~isfield(stimulus_out,'y_start'), stimulus_out.y_start = default.y_start; end

        if ~isfield(stimulus_out,'y_end'), stimulus_out.y_end = default.y_end; end

        if ~isfield(stimulus_out,'monitor_x'), stimulus_out.monitor_x = params.monitor_x; end

        if ~isfield(stimulus_out,'monitor_y'), stimulus_out.monitor_y = params.monitor_y; end

        if ~isfield(stimulus_out,'monitor_refresh'), stimulus_out.monitor_refresh = params.monitor_refresh; end

    
    
    case 'oled'
        
       default = struct;

       % square
        default.x_start = 100;
        default.x_end = 700;
        default.y_start = 0;
        default.y_end = 600;

        default.monitor_refresh = 66.35;
        default.monitor_x = 800;
        default.monitor_y = 600;
        
       % enter missing parameters using the defaults, being careful not to replace any existing fields

        if ~isfield(stimulus_out,'stixel_height')
            stimulus_out.stixel_height = (default.y_end - default.y_start)/stimulus.field_height;
        end

        if ~isfield(stimulus_out,'stixel_width')
            stimulus_out.stixel_width = (default.x_end - default.x_start)/stimulus.field_width;
        end

        if ~isfield(stimulus_out,'x_start'),stimulus_out.x_start = default.x_start; end

        if ~isfield(stimulus_out,'x_end'), stimulus_out.x_end = default.x_end; end

        if ~isfield(stimulus_out,'y_start'), stimulus_out.y_start = default.y_start; end

        if ~isfield(stimulus_out,'y_end'), stimulus_out.y_end = default.y_end; end

        if ~isfield(stimulus_out,'monitor_x'), stimulus_out.monitor_x = default.monitor_x; end

        if ~isfield(stimulus_out,'monitor_y'), stimulus_out.monitor_y = default.monitor_y; end

        if ~isfield(stimulus_out,'monitor_refresh'), stimulus_out.monitor_refresh = default.monitor_refresh; end
        
        
        
end







% show what happened

fprintf('\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
fprintf('!!!  WARNING!  STIMULUS PARAMETERS WERE GUESSED!  !!!\n\n')
fprintf('given parameters:\n\n')
disp(stimulus)
fprintf('guessed parameters:\n\n')
disp(stimulus_out)
fprintf('!!!  WARNING!  STIMULUS PARAMETERS WERE GUESSED!  !!!\n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n')


