function [spatial_sensitivity,all_sig_stixels,cell_ids] = compute_spatial_sensitivity(datarun, cell_spec, varargin)
% compute_spatial_sensitivity     identify stixels with "significant" spatial sampling
%
% usage:  [spatial_sensitivity,all_sig_stixels,cell_ids] = ...
%               compute_spatial_sensitivity(datarun, cell_spec, params)
%
% arguments:  datarun - datarun struct
%           cell_spec - which RFs to use
%              params - struct of optional parameters (see below)
%
% outputs:
%          spatial_sensitivity - Y x X double matrix of spatial sensitivity at each stixel
%              all_sig_stixels - (X*Y) x N boolean matrix of which cells contributed to which stixels
%                                   each column is a cell, each row is a stixel
%                     cell_ids - Nx1 vector, cell ids of the RFs that were used
%
%
% optional fields in params, their default values, and what they specify:
%
% roi                   []          only cells within the ROI are analyzed
% combine_stixels       'max'       how to combine, at each stixel, the sensitivity values of different RFs
%                                       'max' - take the maximum 
%                                       'sum' - take the sum
% foa_spat_sens         []          figure or axes to plot spatial sensitivity as it fills in.
%                                       if 0, make new.  if empty, don't plot
% pause                 []          Seconds to pause after adding each cell
%                                       to plot.  For you to look at it!
% verbose               false       show what's going on
% colors                'bw'        combine RGB?
%                                       'bw' - no
%                                       'rgb' - yes
%
%
% parameters passed on to 'find_cones_in_rf'. if not specified by the user, these parameters are not passed.
%
%       strength                how to summarize RGB values in the RF
%       filter                  how to filter RF before looking for cones
%       selection_params      	struct of params to pass to function find_cones_in_rf
%       fig_single_cell         (passed as fig_rfs)   figure number to plot each RF.
%
%
%
% gauthier 2008-09
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('roi', []);
p.addParamValue('combine_stixels', 'max');
p.addParamValue('foa_spat_sens', []);
p.addParamValue('pause', []);
p.addParamValue('verbose', false);
p.addParamValue('colors', 'bw', @(x)any(strcmpi(x,{'bw','rgb'})));

% parameters to be passed on to find_cones_in_rf
p.addParamValue('strength', 'default value');
p.addParamValue('filter', 'default value');
p.addParamValue('selection_params', 'default value');
p.addParamValue('fig_single_cell', 1);
p.addParamValue('fig_rfs', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% make params struct to pass to find_cones_in_rf
find_params = make_struct_to_pass(p.Results,{'selection','selection_params',...
    'filter','filter','fig_rfs','fig_single_cell','strength','strength'});
find_params.fig_rfs = params.fig_rfs;


% BODY OF THE FUNCTION

% get list of cell numbers
if isempty(params.roi)
    cell_nums = get_cell_indices(datarun, cell_spec);
else
    cell_nums = get_cell_indices_roi(datarun, cell_spec, params.roi);
end

    

% initialize variable to store union of all cones
spatial_sensitivity = [];

% initialize variable that stores significant stixels
all_sig_stixels = false(datarun.stimulus.field_height,datarun.stimulus.field_width,length(cell_nums));


if params.verbose
    start_time = clock; % note when it started
end
    
% set up plot
pa_spat_sens = set_up_fig_or_axes(params.foa_spat_sens);


% Print progress bar
csswb = css_waitbar(params, sprintf('Looking for cones in %d cells...',length(cell_nums)));

% Run analysis
for cc = 1:length(cell_nums)

    % get rf frame
    rf = get_rf(datarun,datarun.cell_ids(cell_nums(cc)));

    % if the rf frame is empty, indicate that there are no significant stixels
    if isempty(rf)
        continue
    end
    
    % get significant stixels for this cell
    [sig_stixels,rf_strength] = find_cones_in_rf(rf,find_params);
    %[sig_stixels,rf_strength] = find_cones_in_rf(rf,struct('fig_rfs',params.fig_single_cell,'selection',params.selection_params));

    % expand
    %sig_stixels = boolean(imfilter(double(sig_stixels),ones(3,3)));
    
    % save sig_stixels
    all_sig_stixels(:,:,cc) = sig_stixels;
    
    % deal with color
    switch params.colors
        case 'bw'
            % do nothing!
        case 'rgb'
            rf_strength = abs(rf);
            sig_stixels = repmat(sig_stixels,[1 1 size(rf,3)]);
    end
    
    % combine this with the existing cone field
    if isempty(spatial_sensitivity)
        % if this is the first cell, make spatial_sensitivity equal to the rf strength
        spatial_sensitivity = sig_stixels.*rf_strength;
    else
        switch params.combine_stixels
            case 'max' % max
                spatial_sensitivity = max(sig_stixels.*rf_strength, spatial_sensitivity);
            case 'sum' % sum
                spatial_sensitivity = sig_stixels.*rf_strength + spatial_sensitivity;
        end
    end
    
    % zero out edges of the cone field
    %border_width = 25;
    %spatial_sensitivity([1:border_width end-border_width+1:end],:) = 0;
    %spatial_sensitivity(:,[1:border_width end-border_width+1:end]) = 0;
    
    % show the growing cone field
    if ~isempty(pa_spat_sens)
        axes(pa_spat_sens)
        switch params.colors
            case 'bw'
                imagesc(spatial_sensitivity)
            case 'rgb'
                imagesc(norm_image(spatial_sensitivity))
        end
        axis image
        drawnow
        
        if ~isempty(params.pause)
            pause(params.pause)
        end
    end
    
    csswb = css_waitbar(params, csswb, cc/length(cell_nums));
end


% reshape all_sig_stixels so that each row is a stixel, and each column is a cell
all_sig_stixels = reshape(all_sig_stixels,[],size(all_sig_stixels,3));



if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

% note which cell ids were used
% this is the same length as the number of columns in all_sig_stixels the matrix all_sig_stixels
cell_ids = datarun.cell_ids(cell_nums);



function csswb = css_waitbar(params, varargin)
csswb = [];
if ~params.verbose, return; end
csswb = text_waitbar(varargin{:});