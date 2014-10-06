function output = fit_cone_rfs(datarun,cell_spec,varargin)
% fit_cone_rfs     fit DOG to single cone RFs
%
%   for each RGC, a fixed center point is chosen, and a circular DOG is fit.  The
%   fitting is eased by first fitting a single gaussian, and basing the DOG on this shape.
%
%
%  TWO MODES:
%
%       1) if cone information is passed in through separate variables,
%           then this function returns just a cell array of fit structs
%
%       2) if cone information is passed in via the datarun struct,
%           then this function returns the datarun struct with the new fits stored in the right place
%
%
%
% usage:  output = fit_cone_rfs(datarun, cell_spec, varargin)
%
% arguments:
%        cone_centers - Nx2 matrix, x,y coordinates of each cone center point
%         rgc_centers - Rx2 matrix, x,y coordinates of each RGC RF center point
%        cone_weights - NxR matrix of cone weights
%            varargin - struct or list of optional parameters (see below)
%
% outputs:
%
%   mode 1:   fits - cell array of struct describing the fit to each RGC
%   mode 2:   datarun - datarun struct with fits stored in datarun.cones.rf_fits
%
%
% optional params, their default values, and what they specify:
%
% center_type       'com'           what kind of center point to use for the fixed center
% fit_radius       	Inf         	radius in which to compute fit
% verbose           false          	show progress ticks
% show_fit_params   false           show parameters of each fit
% foa_profile       []             	figure or axes to plot watch each profile being fit.
%                                       if 0, make new figure. if empty, don't plot
% foa_2d            []             	figure or axes to plot reconstructed RF and the final fit.
%
%
% cone_info         []              if empty, get all cone info from datarun
%                                       otherwise, cone_info should have the following fields:
%
%  	cone_weights	CxR matrix of cone weights (C = #cones, R = #RGCs).  if empty, use datarun.cones.weights
% 	cone_centers	CX2 matrix of center points for each cone.  if empty, use datarun.cones.centers
% 	cell_ids        list of cell ids of the rgcs in cone_weights
%
%
%
% examples: 
%
%   mode 1:
%
%      	cone_info.cone_weights = cone_weights;
%       cone_info.cone_centers = cone_centers;
%       cone_info.cell_ids = cell_ids;
%       rf_cone_fits = fit_cone_rfs(datarun,{2},'cone_info',cone_info);
%
%
%   mode 2:
%
%       datarun = fit_cone_rfs(datarun,{2});
%
%
% gauthier 2008-10
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('center_type', 'com');
p.addParamValue('fit_radius', Inf);
p.addParamValue('verbose', false);
p.addParamValue('show_fit_params', false);
p.addParamValue('foa_profile', []);
p.addParamValue('foa_2d', []);
p.addParamValue('cone_info', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
pa_2d = set_up_fig_or_axes(params.foa_2d);




% get cell indices, cone center points, and initialize variable to store the fits


if isempty(params.cone_info)
    % mode 2
    % if no info was specified, get everything from datarun
    
    cell_indices = get_cell_indices(datarun,cell_spec);
    cone_centers = datarun.cones.centers;
    
    % initialize list fits to be the fits currently stored in datarun
    if isfield(datarun,'cones') && isfield(datarun.cones,'rf_fits')
        fits = datarun.cones.rf_fits;
    else
        % if no fits are stored in datarun, start with an empty cell array
        fits = cell(length(datarun.cell_ids),1);
    end
    
    
else
    % mode 1
    % use cone_info
    
    % read in cone center points
    cone_centers = params.cone_info.cone_centers;
    % get cell_indices, i.e. which cells (columns) to plot in cone_weights
    % this command identifies which columns in cone_weights correspond to the desired cells of cell_spec
    [junk,junk2,cell_indices] = intersect(...
        datarun.cell_ids(get_cell_indices(datarun,cell_spec)), params.cone_info.cell_ids );
    % initialize cell array to store fits
    fits = cell(length(cell_indices),1);
end






% show output
if params.verbose
    fprintf('\nComputing fits for %d RGCs',length(cell_indices));
    start_time = clock; % note when it started
end



% go through each RGC
for cc = 1:length(cell_indices)

    % add tick mark
    if params.verbose,fprintf('.'),end

    
    
    % get cone weights and center point for this RGC
    
    if isempty(params.cone_info)
        % mode 2, use datarun
        cell_index = cell_indices(cc);
        fit_index = cell_index; % where to store the fit
        the_weights = datarun.cones.weights(:,cell_index);
        rf_ctr = rf_center(datarun,datarun.cell_ids(cell_index),params.center_type);
    else
        % mode 1, use cone_info
        cell_index = cell_indices(cc);
        fit_index = cc; % where to store the fit
        the_weights = params.cone_info.cone_weights(:,cell_index);
        rf_ctr = rf_center(datarun,params.cone_info.cell_ids(cell_index),params.center_type);

    end


    % if there's no center point, no fit is computed
    if isempty(rf_ctr)
        fits{cc} = [];
        continue
    end


    % plot reconstructed RF
    if ~isempty(pa_2d)
        % get reconstruction
        cone_rf = cone_rf_reconstructed([datarun.stimulus.field_height datarun.stimulus.field_width],...
            the_weights,cone_centers);
        % plot it
        axes(pa_2d);cla;
        image(norm_image(cone_rf)); axis image; drawnow
        % add center point
        hold on; plot(rf_ctr(1),rf_ctr(2),'r')
    end
    
    

    % get cone weights as a function of distance from the RF center
    [x,y] = rf_cone_profile(the_weights, cone_centers, rf_ctr,'radius',params.fit_radius);

    % if no cones were found, no fit is computed
    if isempty(y)
        fits{cc} = [];
        continue
    end
    
    % put in units of SNR
    noise_sigma = robust_std(the_weights);
    if noise_sigma >0
        y = y/noise_sigma;
    end

    % fit single gaussian to get parameters approximately correct
    fit_params_temp = fit_profile(x,y,'fig_or_axes',[],'verbose',0,'fit_center',false,...
        'fit_surround_scale',false,'surround_scale',0,'center_radius',9);

    % fit a DOG based on these parameters
    [fit_params,y_fit] = fit_profile(x,y,'fig_or_axes',params.foa_profile,'verbose',0,'fit_center',false,...
        'center_radius',fit_params_temp.center_radius,'center_scale',1.4*fit_params_temp.center_scale,...
        'surround_radius',2*fit_params_temp.center_radius,'surround_scale',0.33*fit_params_temp.center_scale,...
        'fit_center_radius',true,'fit_center_scale',true,...
        'fit_surround_radius',true,'fit_surround_scale',true);

    % note fit radius
    fit_params.fit_radius = params.fit_radius;

    % note error
    fit_params.error = norm(y_fit(:,2)-y);
    
    % note which center point was used
    fit_params.center = rf_ctr;
    
    if params.show_fit_params
        disp(fit_params)
    end
    
    
    % plot fit
    if ~isempty(pa_2d)
        % generate points to plot a circle at one center sigma
        temp = [0:100]/100*2*pi;
        circ_pts_x = cos(temp);
        circ_pts_y = sin(temp);
        % get center sigma
        fit_rad = fit_params.center_radius;
        % plot circle at 1 center sigma
        axes(pa_2d);
        plot(rf_ctr(1)+fit_rad*circ_pts_x,rf_ctr(2)+fit_rad*circ_pts_y,'r')
        % pause to behold
        pause(1)
    end

    
    
    % save struct of parameters in the growing cell array of fits
    fits{fit_index} = fit_params;

end


% display how long it took
if params.verbose
    fprintf('\n      done (%0.1f seconds)\n',etime(clock,start_time));
end



% mode 2
if isempty(params.cone_info)

    % store fits in datarun
    datarun.cones.rf_fits = fits;

    % return datarun
    output = datarun;

else
    % mode 1

    % return just the cell array of fits
    output = fits;

end


