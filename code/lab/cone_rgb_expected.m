function rgb_weights = cone_rgb_expected(datarun,params)
% expected_cone_weights     Return expected cone sensitivity to the RGB monitor guns
%
% usage:  rgb_weights = expected_cone_weights(settings,params)
%
% arguments:  datarun - datarun struct with fields
%                           datarun.piece.rig
%                           datarun.piece.optical_path_direction
%                       these are used to identify which rig and optical path was used
%
%              params - struct of optional parameters (see below)
%
%
% outputs:
%         rgb_weights - struct with the following fields
%
%               rgb_weights.L - 1x3 vector of relative expected sensitivity
%                                   to the R, G, and B guns for L cones
%               rgb_weights.M - same for M cones
%               rgb_weights.S - same for S cones
%
% optional fields in params, their default values, and what they specify:
%
% normalize         true        normalize amplitude so that the sum of squares is 1
% figure            []          which figure to plot in.  if 0, make new.  if empty, don't plot.
%
%
%
% gauthier 2008-09
%   NOTE: DATE SPECIFICATION NOT IMPLEMENTED YET
%




% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.normalize = true;
defaults.figure = [];

% combine user and default parameters
params = default_params(defaults, params);


if strcmp(datarun.stimulus.independent, 'nil')
    rgb_weights.C = [1 1 1];
else

    % get rig and path
    rig_path_display = [datarun.piece.rig '-' datarun.piece.optical_path_direction '-' datarun.piece.display];

    % return appropriate weights
    switch rig_path_display

        case 'B-below-crt2'
            rgb_weights.L = [7.1139 15.575	4.0318];
            rgb_weights.M = [2.56	15.814	6.4323];
            rgb_weights.S = [0.1829	0.83108	10.461];

        case 'B-above-crt2'
            rgb_weights.L = [822.97 2533.5 746.26];
            rgb_weights.M = [325.85 2798.8 1213.6];
            rgb_weights.S = [36.843 177.56 2580.8];
            
        case 'B-below-oled2'    
            rgb_weights.L = [2.8807 5.0302 2.0867];
            rgb_weights.M = [1.4706 5.3799 3.0780];
            rgb_weights.S = [0.2943 0.5480 2.1230];
            
        case 'A-below-crt1'
            % !!!!!!!!!  FAKE NUMBERS!!!!!!!!!!!!      !!!!!!!!!!!FAKE NUMBERS!!!!!!!!!!!!!!!!!!!!!!!
            %rgb_weights.L = [0.4294    0.8737    0.2134];
            %rgb_weights.M = [0.1592    0.9193    0.3502];
            %rgb_weights.S = [0.0174    0.1027    0.9618];
            % these are reasonable placeholders, but are not based on measurements
            rgb_weights.L = [20.8575   66.9626   16.4906];
            rgb_weights.M = [ 8.1011   71.6135   25.1743];
            rgb_weights.S = [ 1.0147    2.2233   18.2843];
            disp('!!!!!!!!!!!  USING INCORRECT CONE RGB VALUES   !!!!!!!!!!!')

        case 'A-above-crt1'
            % !!!!!!!!!  FAKE NUMBERS!!!!!!!!!!!!      !!!!!!!!!!!FAKE NUMBERS!!!!!!!!!!!!!!!!!!!!!!!
            % these are not even reasonable placeholders!
            rgb_weights.L = [20.8575   66.9626   16.4906];
            rgb_weights.M = [ 8.1011   71.6135   25.1743];
            rgb_weights.S = [ 1.0147    2.2233   18.2843];
            disp('!!!!!!!!!!!  USING INCORRECT CONE RGB VALUES   !!!!!!!!!!!')

        case 'A-below-oled1'    
            rgb_weights.L = [2.8807 5.0302 2.0867];
            rgb_weights.M = [1.4706 5.3799 3.0780];
            rgb_weights.S = [0.2943 0.5480 2.1230];

        otherwise
            error('expected cone weights for light path ''%s'' not yet computed', rig_path_display)
    end

    % normalize amplitudes, if desired
    if params.normalize
        rgb_weights.L = rgb_weights.L/sqrt(sum(rgb_weights.L.^2));
        rgb_weights.M = rgb_weights.M/sqrt(sum(rgb_weights.M.^2));
        rgb_weights.S = rgb_weights.S/sqrt(sum(rgb_weights.S.^2));
    end


    % plot
    if ~isempty(params.figure)
        % get figure number
        if params.figure == 0
            params.figure = figure;
        end
        % clear it
        figure(params.figure);clf

        % plot bars showing color sensitvity

        subplot(311);bar(rgb_weights.L)
        title('L cone')
        set(gca,'XTickLabel',{'R','G','B'})

        subplot(312);bar(rgb_weights.M)
        title('M cone')
        set(gca,'XTickLabel',{'R','G','B'})
        ylabel('relative sensitivity')

        subplot(313);bar(rgb_weights.S)
        title('S cone')
        set(gca,'XTickLabel',{'R','G','B'})
        xlabel('monitor gun')

        % set size nicely
        curr_size = get(params.figure,'Position');
        set(params.figure,'Position',[curr_size(1) curr_size(2) 230 590])

        drawnow
    end
end


