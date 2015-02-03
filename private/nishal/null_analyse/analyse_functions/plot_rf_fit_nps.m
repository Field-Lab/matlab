function X=plot_rf_fit_nps(datarun, cell_specification, varargin)
%
% plot_rf_fit     plots outlines of fits to spatial receptive field
%
% usage:  h = plot_rf_fit(datarun, cell_specification, varagin)
%
% arguments:     datarun - datarun structure
%     cell_specification - specifies cells of interest
%                          see get_cell_indices for information on ways to specify cells and cell types
%               varargin - struct or list of optional parameters (see below)
%
% outputs:             h - figure handle
%
%
% optional params, their default values, and what they specify:
%
% fits_to_use       'vision'            'vision' or 'obvius' fits
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% sd_radius         1                   contour level in units of sigmal
% edge_color        'k'                 specifies color by letter or color triplet vector (i.e. [0 0 0])                                          
% fill              false               fill the ellipse
% fill_color        'w'                 color of ellipse fill                                         
%
%
% greschner modified from PLOT_RF_OUTLINES shlens 2006-07-25
% gdf: documented and expanded to include obvius fits
% greschner: fill color as matix [index:RGB]


% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addRequired('datarun', @isstruct);
p.addRequired('cell_specification');

% SET UP OPTIONAL ARGUMENTS
p.addParamValue('fits_to_use', 'vision', @ischar);
p.addParamValue('fig_or_axes', -1, @isnumeric);
p.addParamValue('clear_axes', false, @islogical); 

p.addParamValue('sd_radius', 1, @isnumeric);
p.addParamValue('edge_color', [0 0 0]); 
p.addParamValue('edge', true, @islogical);
p.addParamValue('fill', false, @islogical);
p.addParamValue('fill_color', [0 0 0]); 
p.addParamValue('alpha', 1);
p.addParamValue('scale', 1); 
p.addParamValue('labels', false, @islogical);
p.addParamValue('array', false, @islogical)

% removed by gdf
%p.addParamValue('reverse_y', 0);%eg plot combine with image 
%p.addParamValue('scale', 1);%  

% parse inputs
p.parse(datarun, cell_specification, varargin{:});
params = p.Results;

% set fits (vision vs. obvius) to use
if ~isfield(datarun, 'default_sta_fits')
    datarun.default_sta_fits = params.fits_to_use;
else
    default_fits = datarun.default_sta_fits;
    if ~strcmp(default_fits, params.fits_to_use)
        warning('The default fit locations set in datarun, do not match specified fits; using specified fits: ', params.fits_to_use)
        datarun.default_sta_fits = params.fits_to_use;
    end
end
 

if params.fill && size(params.fill_color,1)==1
   params.fill_color=ones(size(datarun.cell_ids))'*params.fill_color;
end



% BODY OF FUNCTION

% get cell numbers
indices = get_cell_indices(datarun,cell_specification);
cell_ids = datarun.cell_ids(indices);

% circle samples
circle_samples = 0:0.005:2*pi;
x_circle = cos(circle_samples);
y_circle = sin(circle_samples);

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes,params.clear_axes);
axes(plot_axes); 
hold on;

% draw each cell
for i=1:length(indices)
  
    % grab parameters
    sd = datarun.(datarun.default_sta_fits).sta_fits{indices(i)}.sd*params.scale;
    mn = datarun.(datarun.default_sta_fits).sta_fits{indices(i)}.mean*params.scale;
    angle = datarun.(datarun.default_sta_fits).sta_fits{indices(i)}.angle;
    
    if strcmp(params.fits_to_use, 'vision')

        %adjust orientation
        angle = -angle;

        % rotate by angle and stretch
        R = rmatrix2(angle / (2*pi) * 360);
        L = params.sd_radius * [sd(1) 0; 0 sd(2)];
        X = R * L * [x_circle; y_circle];
        X(:,end+1) = X(:,1);

        X(1,:)=(X(1,:)+mn(1));
        X(2,:)=(abs((X(2,:)+mn(2))));
    
    elseif strcmp(params.fits_to_use, 'obvius')
    
        % adjust_mean
        mn(2) = datarun.stimulus.field_height - mn(2);

        %adjust orientation
        angle = (pi/2) - angle;

        % rotate by angle and stretch
        R = rmatrix2(angle / (2*pi) * 360);
        L = params.sd_radius * [sd(2) 0; 0 sd(1)];
        X = R * L * [x_circle; y_circle];
        X(:,end+1) = X(:,1);

        X(1,:)=(X(1,:)+mn(1));
        X(2,:)=(abs((X(2,:)+mn(2))));

    else
        error('Unsupported fit specification: must be vision or obvius')
    end
  

    if params.fill %& ~all(params.fill_color(indices(i),:))
        p = patch(X(1,:), X(2,:),  params.fill_color(indices(i),:),'FaceAlpha',params.alpha,'EdgeColor','none');
    end
    if params.edge
        p = plot(X(1,:), X(2,:), 'color', params.edge_color);
    end
    if params.labels
        p = text(mn(1),mn(2), num2str(cell_ids(i)), 'color', params.edge_color);
    end
    
end

% plot array outline
if params.array
    %switch params.coordinates
    %    case 'monitor'
    %        plot_array(datarun)
    %    case 'sta'
            T = coordinate_transform(datarun,'monitor');
            plot_array(datarun,'k',fliptform(T));
    %    otherwise
    %        error('coordinate system ''%s'' not recognized',params.coordinates)
    %end
        
end

%hold off





