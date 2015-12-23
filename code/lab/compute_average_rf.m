function [average_rf,xdata,ydata] = compute_average_rf(datarun, cell_spec, varargin)
% compute_average_rf     compute the average RF.  RFs are registered by supersampling.
%
% usage:  [average_rf,xdata,ydata] = compute_average_rf(datarun, cell_spec, <params>)
%
% arguments:  datarun -
%           cell_spec - 
%            <params> - struct or list of optional parameters (see below)
%
% outputs:       average_rf - 
%             x_data,y_data - coordinates of average_rf in stixel space
%
%
% optional params, their default values, and what they specify:
%
% scale             8                   how much to supersample the RFs
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure.
%                                           if empty, don't plot.  if -1, plot in current.
% avg_size          [15 15]             size of the average RF (in stixels, actual output will be
%                                           larger by a factor of params.scale)
% normalize         'var'               how to normalize the RF before adding to the average
%                                           'var' - set variance of all stixels to 1
%                                           'none' - un-normalized
% verbose           false               show output
%
%
% 2010-05 (!) gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('scale', 8);
p.addParamValue('fig_or_axes', []);
p.addParamValue('avg_size',[1 1] * 15);
p.addParamValue('verbose',false);
p.addParamValue('normalize','var',@(x)any(strcmpi(x,{'var','none'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;
avg_size = params.avg_size;
scale = params.scale;



% initialize 
average_rf = [];
cell_count = 0;


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% add each cell to it
cell_ids = get_cell_ids(datarun,cell_spec);

% show output
if params.verbose
    fprintf('\nAveraging %d cells...',length(cell_ids));
    start_time = clock; % note when it started
end


for cell_id = cell_ids
    
    
    % get the RF
    rf = get_rf(datarun,cell_id);
    center = rf_center(datarun,cell_id);
    
    % skip if problem
    if isempty(rf) || isempty(center) || any(isnan(center)) || ...
            any(center<0) || any(center > size(rf(:,:,1)))
        continue
    end
    
    % ensure center isn't close to an integer without actually being an integer (dont' ask)
    center = round(center*1000)/1000;
    
    % also skip if too close to the edge
    if any((center-avg_size/2) < [1 1]) || any((center+avg_size/2) > fliplr(size(rf(:,:,1))))
        continue
    end
    
    %figure(2);clf;imagesc(norm_image(rf));axis image; hold on;plot(center(1),center(2),'or')
    
    % take subset in the desired region
    roi_x = round(center(1) + [-1 1]*avg_size(2)/2);
    roi_y = round(center(2) + [-1 1]*avg_size(1)/2);
    rf_subset = rf(roi_y(1):roi_y(2),roi_x(1):roi_x(2),:);
    
    % normalize the RF
    switch params.normalize
        case 'var'
            rf_subset = rf_subset/std(reshape(rf_subset,[],1));
        case 'none'
        otherwise
            error('normalization ''%s'' not recognized.',params.normalize)
    end
    
    % interpolate by a factor of scale
    rf_interp = matrix_scaled_up(rf_subset,scale);
    
    % compute new center point
    new_center = center - [roi_x(1) roi_y(1)] + 1;
    new_center = (new_center - 1) * scale + (scale+1)/2;

    % plot pre and post transformation to verify transformation worked as expected
    if 0
        figure(2);clf;
        subplot(121);imagesc(norm_image(rf));axis image; hold on;plot(center(1),center(2),'or')
        subplot(122);imagesc(norm_image(rf_interp));axis image; hold on;plot(new_center(1),new_center(2),'or')
        pause
    end

    % add to the average matrix
    roi_x = round(new_center(1) + [-1 1]*avg_size(2)*scale/2);
    roi_y = round(new_center(2) + [-1 1]*avg_size(1)*scale/2);
    new_rf = rf_interp(roi_y(1):roi_y(2),roi_x(1):roi_x(2),:);
    
    if isempty(average_rf)
        average_rf = new_rf;
    else 
        average_rf = average_rf + new_rf;
    end
    
    % increment cell count
    cell_count = cell_count + 1;

    % plot it
    if ~isempty(plot_axes)
        axes(plot_axes)
        cla
        imagesc(norm_image(average_rf));axis image
        drawnow
    end   
end

% divide by number of cells
if cell_count
    average_rf = average_rf / cell_count;
end

% note x, y data
%T = @(x)(1 + (x - (scale+1)/2) / scale);
T = @(x)x/scale;
xdata = T([-1 1]*avg_size(2)*scale/2);
ydata = T([-1 1]*avg_size(1)*scale/2);


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

