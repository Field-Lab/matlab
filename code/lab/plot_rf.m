function plot_axes = plot_rf(datarun, cell_id, varargin)
% plot_rf_norm     plot a normalized image of a cell's RF, including various bits of info in the title
%
% usage:  plot_rf(datarun, cell_id, <params>)
%
% arguments:   datarun - argument 1
%              cell_id - which cell
%             <params> - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% foa           0      	figure or axes to plot in. if 0, make new figure. if empty, don't plot
% coordinates   'sta'   what coordinates to plot in
% scale         1       how much to scale up the matrix
% com           false    if available in datarun.stas, plot the center of mass
% sig_stix      false    if available in datarun.stas, plot the significant stixels
% polarity      false   if true and polarity is available, then it is used to plot rf
% array         false   if true, plot array outline
% fit           true    plot rf fit, if available
% color_transform []    change color space
% autozoom      false   Autozoom in based on fit bounds and AZ_PAD_FACTOR
% az_pad_factor 2.5     The factor to scale up the fit bounds for autozoom
% title         true    Show title?
% pdtitle       true    Include piece and datarun info in title
%
% 2009-04 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', -1);
p.addParamValue('coordinates', 'sta', @(x) any(strcmpi(x,{'sta','monitor'})));
p.addParamValue('scale', 1);
p.addParamValue('com', false, @islogical);
p.addParamValue('sig_stix', false, @islogical);
p.addParamValue('polarity', true, @islogical);
p.addParamValue('color_transform', [], @isnumeric)
p.addParamValue('array', false, @islogical);
p.addParamValue('fit', true, @islogical);
p.addParamValue('autozoom', false);
p.addParamValue('az_pad_factor', 5);
p.addParamValue('az_aspect_ratio', 1);
p.addParamValue('title', true);
p.addParamValue('pdtitle', true);
p.addParamValue('ticks', false);
p.addParamValue('gun', [])

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% BODY OF THE FUNCTION


% generate title text
title_text = sprintf('cell id %d',cell_id);

% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa);

% get the RF
rf = get_rf(datarun, cell_id, 'color_transform', params.color_transform);

% scale up
rf = matrix_scaled_up(rf,params.scale);

% identify coordinates to plot in
switch params.coordinates
    case 'sta'
        %xx = [1 size(rf,2)];
        %yy = [1 size(rf,1)];
        tform = coordinate_transform(datarun,'sta','input','sta scaled','scale',params.scale);
    case 'monitor'
        tform = coordinate_transform(datarun,'monitor','input','sta scaled','scale',params.scale);
end
[xx,yy] = tformfwd(tform,[1 size(rf,2)],[1 size(rf,1)]);

% get polarity
if params.polarity
    cell_index  = get_cell_indices(datarun, cell_id);
    if isfield(datarun.stas, 'polarities') && ~isempty(datarun.stas.polarities{cell_index})
        polarity = datarun.stas.polarities{cell_index};
        
        if polarity == 0
            polarity = 1;
        end
        
        % correct RF polarity before plotting
        rf = rf * polarity;
   
        % note change in title text
        title_text = [title_text ' (polarity corrected)'];
    else
        % otherwise, warn that polarity information was missing
        warning('polarity information unavailable, OFF cells might look like ON cells')
    end
end


% plot the image
if ~isempty(params.gun) % check to see if a single primary (gun) is being called for plotting
    imagesc(xx,yy,norm_image(rf(:,:,gun)),'parent',plot_axes);
    axis image; hold on   
else
    imagesc(xx,yy,norm_image(rf),'parent',plot_axes);
    axis image; hold on    
end

% plot fit, if possible
if params.fit
    if isfield(datarun,'stas') && isfield(datarun.stas,'fits')
        plot_rf_summaries(datarun,cell_id,'clear',0,'foa',plot_axes,'plot_fits',1,'coordinates',params.coordinates)
    end
end  


% add cell type, if possible
ct = find_cell_types(datarun,cell_id);
if ct; title_text = [title_text ', ' datarun.cell_types{ct}.name];  end

% plot COM, if possible
if params.com
    ctr = rf_center(datarun,cell_id,'com',params.coordinates);
    if ~isempty(ctr) || ~all(ctr==-1)
        plot(ctr(1),ctr(2),'.r','MarkerSize',10);
        title_text = [title_text ', ' sprintf('COM: (%0.1f, %0.1f)',ctr(1),ctr(2))];
    end
end

% plot significant stixels, if given in datarun
if params.sig_stix
    if isfield(datarun.stas, 'marks')
        cell_index = get_cell_indices(datarun, cell_id);
        %sig_stix = datarun.stas.significant_stixels{cell_index};
        sig_stix = datarun.stas.marks{cell_index};
        temp_indices = find(full(sig_stix));
        [temp_rows, temp_cols] = ind2sub(size(sig_stix), temp_indices);
        plot(temp_cols, temp_rows, 'ko')
    end
end

% plot array outline
if params.array
    switch params.coordinates
        case 'monitor'
            plot_array(datarun)
        case 'sta'
            T = coordinate_transform(datarun,'monitor');
            plot_array(datarun,'k',fliptform(T));
        otherwise
            error('coordinate system ''%s'' not recognized',params.coordinates)
    end
        
end

% Autozoom?
if params.autozoom && isfield(datarun,'stas') && isfield(datarun.stas,'fits')
    autozoom_to_fit(datarun, cell_id, params.az_pad_factor, 1, params.az_aspect_ratio);
end

    
if params.pdtitle
    title_text = {strrep(datarun.names.short_name, '_', '\_'); title_text};
end

% make title
if params.title
    title(title_text)
end

if ~params.ticks
    set(plot_axes, 'XTick', [], 'YTick', []);
end