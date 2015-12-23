function h = plot_cell_sampling(datarun, cell_spec, varargin)
% plot_cell_sampling     make "spider" plot or contour plot
%
%   cells plotted are the interection of those represented in cone_weights and those specified in cell_spec
%
%
% usage:  h = plot_cell_sampling(datarun, cell_spec, params)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells to plot (NOTE: only cells represented in cone_weights will be plotted)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:          h - handle to the plot axes
%
%
% optional fields in params, their default values, and what they specify:
%
% fig_or_axes       0               which figure or axes to plot in
% clear             true            clear axes before proceeding
% plot_cones        true            plot cone mosaic
% roi               []              ROI in which to plot cells.  if empty, plot all.
% label             false           plot the names of cell ids
% fits              false           plot single cone fits
% pause             false           pause after each cell is plotted
%
% cell_colors       []              option to plot each cell a different color
%                                       [] - plot all cell in black
%                                       Nx3 matrix - plot cells according to the specified colors
%                                           if there are more than N cells, the colors wrap back to the beginning
%                                       'net_input' - plot color of each cell to reflect its net L and M input.
%                                          	more opponency means brighter color
%                                       'mix' - plot each cell a different random color
%
% type              'spider'        how to plot each cell
%                                       'spider' - line connecting cell center to each cone
%                                      'contour' - patch describing contour level of RF with cones blurred
%
%   additional options for plot_cones == 'spider'
%       line_width          [realmin 0.6]
%                                       width of the lines
%       center_type         'com'       what type of center point
%       negative_scale      1           scale the line width of negative weights with this factor
%
%   additional options for plot_cones == 'contour'
%       alpha               1          	alpha of contour patches
%
%
%
%
%
%
%
% these parameters are only passed to plot_cone_mosaic if specified by the user:
%
%       cone_size           size of each cone's dot
%     	bg_color            color of the background
%
%
% these parameters are only passed to get_net_cone_input if specified by the user:
%
%       input_radius        radius in which to get the net cone input
%     	radius_type         units of the radius
%     	normalize           how to normalize cone weights
%
%
% these parameters are only passed to select_cone_weights if specified by the user:
%
%       thresh
%       plot_radius
%       polarity            by default, 1 is passed to only plot positive weights
%       contiguity
%       scale
%
%
% gauthier   2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fig_or_axes', 0);
p.addParamValue('clear', true);
p.addParamValue('plot_cones', true);
p.addParamValue('cell_colors', []);
p.addParamValue('roi', []);

p.addParamValue('type', 'spider',@(x)any(strcmp(x,{'contour','spider'})));

% add other parameters based on 'type'
p.KeepUnmatched = true;
p.parse(varargin{:});
switch p.Results.type
    case 'spider'
        p.addParamValue('line_width', [0.1 0.6]);
        p.addParamValue('center_type', 'com');
        p.addParamValue('negative_scale', 1);
    case 'contour'
        p.addParamValue('alpha', 1,@(x)(x>0)|(x<1));
end
p.KeepUnmatched = false;


p.addParamValue('label', false);
p.addParamValue('fits', false);
p.addParamValue('pause', false);


% parameters to be passed on plot_cone_mosaic
p.addParamValue('cone_size', 'default value');
p.addParamValue('bg_color', 'default value');

% parameters to be passed on get_net_cone_input
p.addParamValue('input_radius', 'default value');
p.addParamValue('radius_type', 'default value');
p.addParamValue('normalize', 'default value');

% parameters to be passed on select_cone_weights
p.addParamValue('thresh', 'default value');
p.addParamValue('polarity', 1);
p.addParamValue('plot_radius', 'default value');
p.addParamValue('contiguity', 'default value');
p.addParamValue('scale', 'default value');



% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
mosaic_params = make_struct_to_pass(p.Results,{'cone_size','cone_size','bg_color','bg_color'});
net_input_params = make_struct_to_pass(p.Results,{'radius','input_radius','radius_type','radius_type','normalize','normalize'});
selection_params = make_struct_to_pass(p.Results,{'thresh','thresh','polarity','polarity','radius','plot_radius',...
                    'contiguity', 'contiguity', 'scale','scale'});



% BODY OF THE FUNCTION


if params.plot_cones
    % plot cone mosaic
    plot_axes = plot_cone_mosaic(datarun, mosaic_params,'clear',params.clear,'fig_or_axes',params.fig_or_axes);
else
    % just set up axes
    plot_axes = set_up_fig_or_axes(params.fig_or_axes,params.clear);
end


% get list of cell indices
cell_indices = get_cell_indices_roi(datarun,cell_spec,params.roi);

% set up cell colors
if isempty(params.cell_colors)
    params.cell_colors = [0 0 0];
else
    if strcmp(params.cell_colors,'net_input')
        % get net cone input of each cell, if desired
        net_inputs = get_net_cone_input(datarun,cell_spec,net_input_params);

        % get opponency index of each cell
        opponency_indices = atan2(net_inputs(:,1),net_inputs(:,2)) - pi/4;

        % determine colors
        params.cell_colors = zeros(length(cell_indices),3);
        % non opp
        %params.cell_colors(1:length(cell_indices),:) = repmat([0 0 0],length(cell_indices),1);
        % L excit, M inhib
        temp=find(opponency_indices>=0);
        %params.cell_colors(temp,:) = repmat([1 0 0],length(temp),1);
        params.cell_colors(temp,:) = [min(opponency_indices(temp)/(3*pi/4),1) zeros(length(temp),2)];
        % M excit, L inhib
        temp=find(opponency_indices<0);
        %params.cell_colors(temp,:) = repmat([0 1 0],length(temp),1);
        params.cell_colors(temp,:) = [zeros(length(temp),1)...
            min(-opponency_indices(temp)/(3*pi/4),1) zeros(length(temp),1)];
    end
    if strcmp(params.cell_colors,'mix')
        % note the current random seed
        current_rand_seed = rand('twister');
        % get mix of colors
        cols = spring(length(cell_indices));
        % generate random colors (that will be fixed)
        rand('state',11111);
        params.cell_colors = cols(randperm(size(cols,1)),:);
        % restore seed
        rand('twister',current_rand_seed);
    end
end



% plot all cells together on cone mosaic


% go to plot axes
axes(plot_axes)
hold(plot_axes,'on')

% go through list of cells
for cc = 1:length(cell_indices)
    
    % index of the cell in cell_ids
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % select cone weights
    [the_weights,selection] = select_cone_weights(datarun,cell_id, selection_params);
 
 
    % if none were selected, skip
    if all(~selection)
        fprintf('cell id %d skipped because no weights were found\n',cell_id)
        continue
    end
    
    % get color for this cell
    cell_color = params.cell_colors(mod(cc-1,size(params.cell_colors,1))+1,:);
    %cell_color = 'k';


    switch params.type

        case 'spider'  % draw lines showing which cones are sampled

            % get cell center point
            rf_ctr = rf_center(datarun,datarun.cell_ids(cell_index),params.center_type);

            % skip if no center point
            if isempty(rf_ctr); 
                fprintf('cell id %d skipped because no center point was found\n',cell_id)
                continue
            end

            % get indices of the weights to be plotted
            wt_ind = find(selection);
            

            % compute the line widths based on the weights

            % interpret params.line_width
            plot_width_range = max(params.line_width)-min(params.line_width);
            plot_width_min = min(params.line_width);

            % identify negative cones
            wt_ind_minus = the_weights(wt_ind) < 0;

            % scale them up
            the_weights(wt_ind(wt_ind_minus)) = the_weights(wt_ind(wt_ind_minus)) * params.negative_scale;

            % make all weights positive
            the_weights = abs(the_weights);

            % normalize weights amplitude
            the_weights = the_weights / max(the_weights(selection));

            % compute line widths
            widths = (the_weights(selection)*plot_width_range) + plot_width_min;


            % plot negative cones
            for ww = find(wt_ind_minus)'

                % get cone center point
                cone_ctr = datarun.cones.centers(wt_ind(ww),:);

                % plot line
                plot([rf_ctr(1) cone_ctr(1)],[rf_ctr(2) cone_ctr(2)],'Color',1-cell_color,'LineWidth',widths(ww))
            end
            
            
            % plot positive cones
            for ww = find(~wt_ind_minus)'

                % get cone center point
                cone_ctr = datarun.cones.centers(wt_ind(ww),:);

                % plot line
                plot([rf_ctr(1) cone_ctr(1)],[rf_ctr(2) cone_ctr(2)],'Color',cell_color,'LineWidth',widths(ww))
            end
            
            
            
            %             % draw line from center point to each cone
            %             for ww = 1:length(wt_ind)
            %
            %                 % get cone center point
            %                 cone_ctr = datarun.cones.centers(wt_ind(ww),:);
            %
            %                 % identify if cone is negative
            %                 if wt_ind_minus(ww)
            %                     cell_color_temp = 1-cell_color;
            %                 else
            %                     cell_color_temp = cell_color;
            %                 end
            %
            %                 % plot line
            %                 plot([rf_ctr(1) cone_ctr(1)],[rf_ctr(2) cone_ctr(2)],'Color',cell_color_temp,'LineWidth',widths(ww))
            %             end

            
        case 'contour'  % plot patch around the selected cones
            
            switch 2
                case 1 % convex hull
                    % skip if few than 3 strong cones (convhull can't handle it)
                    if sum(selection)<=2;continue;end

                    % get center points of selected weights
                    x = datarun.cones.centers(selection,1);
                    y = datarun.cones.centers(selection,2);

                    % compute convex hull
                    K = convhull(x,y);
                    hull = [x(K) y(K)];

                    % plot it
                    plot(gca,hull(:,1),hull(:,2),'-','Color',cell_color,'LineWidth',0.5)
                    
                case 2 % contours
                    
                    % get reconstructed RF
                    rf = cone_rf_reconstructed([datarun.stimulus.field_height datarun.stimulus.field_width],...
                        the_weights,datarun.cones.centers,'cones',find(selection));
                    
                    % skip if all zeros
                    if all(all(all(rf==0))); continue; end
                    
                    switch 2
                        case 1 % threshold for all kept cones
                            rf = double(rf~=0);
                            
                            filt = make_gaussian('x_size',41,'y_size',41,'center_radius',1,'center',[21 21]);
                            
                        case 2 % don't threshold
                            
                            %filt = make_gaussian('x_size',41,'y_size',41,'center_radius',3,'center',[21 21]);
                            filt = make_gaussian('x_size',41,'y_size',41,'center_radius',1.5,'center',[21 21]);
                    end
                    
                    % blur
                    rf_blurred = imfilter(rf,filt);
                    
                    % normalize
                    rf_norm = rf_normalized(rf_blurred);
                    
                    % get contours
                    contour_polygons = rf_contours(rf_norm, [2]);
                    
                    
                    % plot contours
                    
                    % set colors
                    cnt_col = {'k'};
                    
                    % go through each threshold
                    for tt = 1:length(contour_polygons)
                        % get color
                        col = cnt_col{mod(tt-1,cc) + 1};

                        % within each threshold, go through each contour band
                        for bb = 1:length(contour_polygons{tt})
                            patch('XData',contour_polygons{tt}(bb).x,'YData',contour_polygons{tt}(bb).y,...
                                'FaceColor',cell_color,'Parent',plot_axes,'LineStyle','none','FaceAlpha',params.alpha);
                            %plot(contour_polygons{tt}(bb).x,contour_polygons{tt}(bb).y,'color',cell_color)
                        end
                    end
                    
                    drawnow
                    
            end

    end

    % pause
    if params.pause
        pause
    end

end


% plot fits
if params.fits
    plot_rf_summaries(datarun,cell_spec,'clear',0,'label',0,'foa',plot_axes,'plot_fits',1)
end
    
% label cells
if params.label
    plot_rf_summaries(datarun,cell_spec,'foa',plot_axes,'label',1,'clear',0,...
        'label_color',[0 0 0],'label_size',15)
    
end



% return handle, if desired
if nargout > 0
    h = plot_axes;
end

