function avg_rf = get_average_rf(datarun, cell_spec, varargin)
% get_average_rf     Average a collection of RFs together.
%
%   before averaging, each RF is linearly interpolated using the provided
%   scale factor, and re-centered.  
% 
%
% usage:  avg_rf = get_average_rf(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     avg_rf - average Rf
%
%
% optional params, their default values, and what they specify:
%
% verbose           false           show text output
% fig_or_axes       []              figure or axes to plot in. if 0, make new figure. if empty, don't plot
% scale             3               linear scale factor for interpolating each RF
% normalize         true            normalize each RF to have a maximum value of 1
% center            []              center point (in stixel coordinates) where to align the RFs
%                                       if empty, uses datarun.stimulus to find the center
% which_rfs         'rfs'           name of field in datarun.stas where the RFs are stored
%
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('scale',3,@(x)(x>=1)&&(mod(x,1)==0));
p.addParamValue('normalize', true);
p.addParamValue('center', []);
p.addParamValue('which_rfs', 'rfs');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% identify location to which all RFs will be averaged
if isempty(params.center)
    new_ctr = [datarun.stimulus.field_width datarun.stimulus.field_height]/2;
else
    new_ctr = params.center;
end

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% show output
if params.verbose
    fprintf('\nComputing average of %d RFs',length(cell_indices));
    start_time = clock; % note when it started
end


% loop through cells
for cc = 1:length(cell_indices)
    
    if params.verbose;fprintf('.');end

    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);

    % get RF
    rf = get_rf(datarun,cell_id,'where',params.which_rfs);
    
    % normalize RF, if desired
    if params.normalize
        rf = rf/max(max(max(rf)));
    end
        
    % get RF center
    rf_ctr = rf_center(datarun,cell_id);

    % skip if empty
    if isempty(rf) || isempty(rf_ctr);continue;end

    % generate RF centered at center point
    new_rf = scale_and_re_center_rf(rf,rf_ctr,new_ctr,params.scale);
    
    % add to growing accumulation
    if exist('avg_rf','var')
        avg_rf = avg_rf + new_rf;
    else
        avg_rf = new_rf;
    end

    % plot, if desired
    if ~isempty(plot_axes)
        image(norm_image(avg_rf),'Parent',plot_axes)
        axis image
        title(sprintf('average of %d RFs',cc))
        drawnow
    end
end




% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end



function new_rf = scale_and_re_center_rf(rf,rf_ctr,new_ctr,scale)


% note width and height of rf
height = size(rf,1);
width = size(rf,2);
    
% identify center points in rescaled coordinates
rf_ctr_scaled = (rf_ctr - 1) * scale + 1;
new_ctr_scaled = (new_ctr - 1) * scale + 1;

% interpolate RF (if needed)
if scale == 1
    rf_scaled = rf;
else
    
    % generate needed variables
    [X,Y] = meshgrid(1:width,1:height);
    [XI,YI] = meshgrid(1:1/scale:width+1,1:1/scale:height+1);
    XI = XI(1:end-1,1:end-1);
    YI = YI(1:end-1,1:end-1);

    % interpolate points, and put a 0 in for any points outside the original range
    for cc = 1:size(rf,3)
        rf_scaled(:,:,cc) = interp2(X,Y,rf(:,:,cc),XI,YI,'linear',0);
    end
end


diff = new_ctr_scaled - rf_ctr_scaled;

xform = [ 1  0  0; 0  1  0; diff  1 ];
tform_translate = maketform('affine',xform);

new_rf = imtransform(rf_scaled, tform_translate,...
    'XData',[1 size(rf_scaled,2)],'YData', [1 size(rf_scaled,1)]);


