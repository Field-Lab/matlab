function qd_plot_aligned_rgc(datarun,axons,cell_id,axon_ids,varargin)
% qd_plot_aligned_rgc     plot an EI and axon, optionally an image in the background
%
% usage:  qd_plot_aligned_rgc(datarun,axons,cell_id,axon_ids,<params>)
%
% arguments:  datarun - argument 1
%               axons - cell array of axon paths, each element an Nx2 matrix of coordinates
%             cell_id - which EI to plot
%            axon_ids - which axon to plot
%            varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% ei_scale          1         	scale on EI circles
% ei_color          [.2 .2 1]   color for ei
% ei_cutoff         -1          cutoff for plotting an electrode's value
% plot_electrodes   false       plot electrodes
% foa               0           figure or axes to plot in. if 0, make new figure. if empty, don't plot.  if -1, plot in current.
% clear             true        clear axes before plotting?
% title             ''          title
% xlim              []          x bounds of plot, if empty use datarun.ei.array_bounds_x
% ylim              []          y bounds of plot, if empty use datarun.ei.array_bounds_y
% bg_color          [1 1 1]     color of background, if empty don't plot background color
%
%
% 2009-08 gauthier
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('ei_scale', 1);
p.addParamValue('ei_color', [.2 .2 1]);
p.addParamValue('ei_cutoff', -1);
p.addParamValue('plot_electrodes', false);
p.addParamValue('foa', 0);
p.addParamValue('clear', true);
p.addParamValue('title', '');
p.addParamValue('xlim',[]);
p.addParamValue('ylim',[]);
p.addParamValue('bg_color', [1 1 1]);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% variable names
ep = datarun.ei.position;
if isempty(params.xlim); params.xlim = datarun.ei.array_bounds_x; end
if isempty(params.ylim); params.ylim = datarun.ei.array_bounds_y; end



% set up plot axes
if params.foa==0
    % if new figure, make subplot
    figure
    plot_axes = subplot('Position',[.05 .02 .95 .98]);
else
    % otherwise, do usual set up
    plot_axes = set_up_fig_or_axes(params.foa,params.clear);
end
axes(plot_axes)

% plot white background the size of the array
if ~isempty(params.bg_color)
    imagesc(permute(params.bg_color,[1 3 2]),'xdata',datarun.ei.array_bounds_x,'ydata',datarun.ei.array_bounds_y)
end

% set up plot parameters
axis image;axis xy;colormap gray;hold on;
set(plot_axes,'ydir','reverse','xdir','reverse','xlim',params.xlim,'ylim',params.ylim)



% plot electrodes
if params.plot_electrodes
    plot(ep(:,1),ep(:,2),'.','Color',[1 1 0])
    for ee=1:size(ep,1)
        text(ep(ee,1),ep(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',10,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
end


% initialize the title text
title_text = params.title;

% plot axon path
for axon_id = axon_ids
    title_text = [title_text 'axon'];
    if ~isempty(axons{axon_id})
        points = traced_cell_points(axons{axon_id}(1:2,:),axons{axon_id}(2:end,:));
        plot(points(:,1),points(:,2),'r')
        % add axon name to title text
        title_text = [title_text sprintf(' %d',axon_id)];
    end
end

% plot EI
if ~isempty(cell_id) && 1
    datarun = load_ei(datarun,cell_id);
    plot_ei(datarun,cell_id,'coordinates','array','alpha',0,'pretty_axes',0,'max_scale',30,...
        'pos_color',params.ei_color,'scale',params.ei_scale,'cutoff',params.ei_cutoff);
    %plot_ei_(new_ei{1},datarun.ei.position,0,'alpha',0,'pretty_axes',0,'pos_color',[1
    %0 0]);
        title_text = [title_text sprintf(', EI %d',cell_id)];
end

% add title
title(title_text);
drawnow();


