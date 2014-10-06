function compare_fits_to_rfs(datarun, varargin)
% COMPARE_FITS_TO_RFS     Examine a 1x1 stixel STA next to the single cone
% fits to that RGC. All cells present in a DATARUN struct will be plotted,
% one at a time. A scrollbar at the bottom of the window is used to
% specify the RGC to examine.
%
% usage:  result = compare_fits_to_rfs(datarun, <params>)
%
% arguments:  datarun - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:       none - the comparison interface is displayed on the screen
%
% optional params, their default values, and what they specify:
%
% figure              []        figure or axes to plot in. if empty make new figure
% cell_specification  'all'     which subset of the cells in datarun should be plotted
% plot_radius         10        plot this pixel radius around each cone
% plot_color          []        plot R,G, and/or B components of images
% cones               true      plot cones on top of images
% cone_size           10        size of cone centers
% scale_factor        1         factor used to scale up STAs/fit images
%
% 4/29/2009 tamachado
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

% specify list of optional parameters
p.addParamValue('cell_specification','all');
p.addParamValue('plot_radius', 10);
p.addParamValue('plot_color', []);
p.addParamValue('figure', []);
p.addParamValue('rfs_field_name', 'rfs');
p.addParamValue('scale_factor', 1);
p.addParamValue('cones', true);
p.addParamValue('cone_size', 10);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the interface
%%%%%%%%%%%%%%%%%%%%%%%%

% set up plot axes
fig = set_up_fig_or_axes(params.figure);

% get cell number
cell_nums = get_cell_indices(datarun,params.cell_specification);

% compute fits
fits = datarun.Wc*datarun.cones.weights;
nCells = size(fits,2);
sy = datarun.stimulus.field_height;
sx = datarun.stimulus.field_width;
fits_reshaped = reshape(fits,[sy sx 3 nCells]);

% set up slider-plotting
k=1; kmin=1; kmax=nCells;

% slider control
ha = uicontrol(gcf,...
    'Style','slider',...
    'Min' ,kmin,'Max',kmax, ...
    'Units','normalized', ...
    'Position',[0,0,.6,.04], ...
    'Value', k,...
    'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
    'CallBack', 'uiresume;');

% toggle button to show/hide cones
conetoggle = uicontrol('Style','togglebutton', ...
   'Units','normalized', ...
   'Position',[0.61 0 .15 .04], ...
   'String','Show Cones',...
   'Callback','uiresume;');
set(conetoggle,'Value',params.cones);

% toggle button to show negative weights only
negtoggle = uicontrol('Style','togglebutton', ...
   'Units','normalized', ...
   'Position',[0.78 0 .15 .04], ...
   'String','Hide Positive',...
   'Callback','uiresume;');
set(negtoggle,'Value',0);

% slider to change plot radius
radspinner = uicontrol(gcf,...
    'Style','slider',...
    'Min' , 1,'Max', max(sx,sy), ...
    'Units','normalized', ...
    'Position',[0,.25,.04,.5], ...
    'Value', k,...
    'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
    'CallBack', 'uiresume;');
set(radspinner,'Value',params.plot_radius);
defaultRadius = params.plot_radius;

while k
    % update gui parameters
    k=round(get(ha,'Value'));
    params.plot_radius  = round(get(radspinner,'Value'));
    params.cones  = get(conetoggle,'Value');
    negativeMode  = get(negtoggle,'Value');
    
    % get the current cell's id
    ind = cell_nums(k);

    % get summary
    summary = get_rf(datarun,datarun.cell_ids(ind));
    fit = fits_reshaped(:,:,:,ind);
    
    % if no color specified, plot all colors
    if ~isempty(params.plot_color)
        % otherwise, plot only the specified color
        summary = summary(:,:,params.plot_color);
        fit = fits_reshaped(:,:,params.plot_color,ind);
    end
    
    % if specified, show only negative weights
    if negativeMode
        fit(fit > 0) = 0;
        summary(summary > 0) = 0;
    end

    % get center point, rounded
    ctr = round(datarun.stas.rf_coms{ind});

    % get plot radius
    rad = params.plot_radius;

    % get plot range, clipping edges if it goes beyond the boundaries of the summary itself
    xrng = max(ctr(1)-rad,1):min(ctr(1)+rad,size(summary,2));
    yrng = max(ctr(2)-rad,1):min(ctr(2)+rad,size(summary,1));

    % plot the STA
    subplot(1,2,1);
    
    % scale up if desired
    if params.scale_factor == 1
        image(norm_image(summary(yrng,xrng,:)))
    else
        image(matrix_scaled_up(norm_image(summary(yrng,xrng,:)),params.scale_factor))
    end
    
    % make a title and rescale the image to be square
    axis image
    type = find_cell_types(datarun,datarun.cell_ids(ind));
    title(sprintf('id: %d, type: %s',datarun.cell_ids(ind), datarun.cell_types{type}.name))
    set(subplot(1,2,1),'XTick',[],'YTick',[])
    
    % plot the fits
    subplot(1,2,2);
    
    % scale up if desired
    if params.scale_factor == 1
        image(norm_image(fit(yrng,xrng,:)));
    else
        image(matrix_scaled_up(norm_image(fit(yrng,xrng,:)),params.scale_factor))
    end
    
    axis image;
    if params.cones,  sc = 'show cones'; else sc = 'hide cones'; end
    if ~negativeMode, af = 'all fits'; else af = 'negative fits only'; end
    title(sprintf('%s, %s',sc, af))
    set(subplot(1,2,2),'XTick',[],'YTick',[])
    
    % plot cones
    if params.cones
        
        % convert cone center points to current coordinates
        cone_centers = datarun.cones.centers;
        cone_centers = cone_centers - repmat([min(xrng)-1 min(yrng)-1],size(cone_centers,1),1);
        cone_centers = params.scale_factor * (cone_centers - 0.5) + 0.5;
        
        % only find cones within radius (so stuff plots super fast)
        ci = intersect(...
             intersect(find(datarun.cones.centers(:,1) >= min(xrng)),...
             find(datarun.cones.centers(:,1) <= max(xrng))),...
             intersect(find(datarun.cones.centers(:,2) >= min(yrng)),...
             find(datarun.cones.centers(:,2) <= max(yrng))));
        

        for sp = 1:2
            subplot(1,2,sp);
            hold on
            for c = ci'
                switch datarun.cones.types(c)
                    case 'L'
                        plotColor = [255 0 102]/255;
                    case 'M'
                        plotColor = [51 255 0]/255;
                    case 'S'
                        plotColor = [0 0 228]/255;
                    case 'U'
                        plotColor = [0 0 0]/255;
                end
                plot_params = struct;
                dotScale = 3; % dots are 3x smaller than the point size
                plot_params.MarkerSize = (defaultRadius/(.5*params.plot_radius))*params.cone_size*dotScale;
                plot(cone_centers(c,1),cone_centers(c,2),'.','Color',plotColor,plot_params)
            end
            hold off
        end
    end

    uiwait;
end