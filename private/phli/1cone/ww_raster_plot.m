function [rax,histx,histy,res] = ww_raster_plot(datarun, stimstruct, cells, polarities, varargin)
% WW_RASTER_PLOT
% usage: [rax,histx,histy,res] = ww_raster_plot(datarun, stimstruct, cells, polarities, varargin)
%
% 2011 phli
%

opts = inputParser();

% These are split out and passed to RASTER
opts.addParamValue('plot', true);
opts.addParamValue('hist', true);
opts.addParamValue('hist_bin', 0.05);
opts.addParamValue('hist_line_width', 2);
opts.addParamValue('color', 0.5*[1 1 1]);
opts.addParamValue('MarkerSize', 5);
opts.addParamValue('FontSize', 10);
opts.addParamValue('start', -0.1);
opts.addParamValue('stop', []);

% Unmatched will also be passed along to raster
opts.KeepUnmatched = true;

% These are mine
opts.addParamValue('triggers', datarun.triggers);
opts.addParamValue('link_ylims', 'all');
opts.addParamValue('set_ylims', true);
opts.addParamValue('hist_colors', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('subplot_add_width', 0);
opts.addParamValue('subplot_add_height', 0);
opts.addParamValue('figure_by_cell_id', false);
opts.addParamValue('plot_repr', false);
opts.addParamValue('titles', false);

opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
if isempty(opts.stop)
    opts.stop = mean(diff(opts.triggers));
end;
rasteropts = keepfields(opts, 'plot', 'hist', 'hist_bin', 'hist_line_width', 'color', 'MarkerSize', 'FontSize', 'start', 'stop');


ncones = stimstruct.numcones;


if length(polarities) == 1
    polarities = repmat(polarities, size(cells));
end


% Get colors
if isempty(opts.hist_colors)
    opts.hist_colors = opts.cmf(ncones);
end


% Get only trials that were actually run (don't get them past number of triggers)
pulses = stimstruct.pulses(1:length(opts.triggers));


% Identify trials for different conditions
intensity = 0.48;
for i = 1:ncones
    up(:,i)    = cell2mat(collect(pulses, @(s) (isequal(s.rgbs{i}, num2cell(intensity * [ 1  1  1])))));
    down(:,i)  = cell2mat(collect(pulses, @(s) (isequal(s.rgbs{i}, num2cell(intensity * [-1 -1 -1])))));
    blank(:,i) = cell2mat(collect(pulses, @(s) (isequal(s.rgbs{i}, num2cell([0 0 0])))));
end
singles = sum(blank,2) == (ncones - 1);
doubles = sum(blank,2) == (ncones - 2);


subplot_width  = ncones + opts.subplot_add_width;
subplot_height = ncones + opts.subplot_add_height;
for i = 1:length(cells)
    cell = cells(i);
    if iscell(cell), cell = cell{1}; end
    if isempty(cell), continue; end
    
    % Set things up according to polarity of cell
    polarity = polarities(i);
    if polarity == 1
        stim = up;
        repr = down;
    elseif polarity == -1
        stim = down;
        repr = up;
    else        
        warning(['Polarity not set for cell ' num2str(cell)]);
        continue
    end
    
    if rasteropts.plot
        if opts.figure_by_cell_id
            figure(cell);
        else
            figure;
        end
    end
    
    % Plot the rasters
    for p = 1:ncones
        for q = 1:ncones
            if rasteropts.plot
                sanesubplot(subplot_height, subplot_width, {p q});
            end
            
            if p == q
                trials = stim(:,p) & singles;
                [res{p,q}, ax, histx, histy{p,q}] = rasterphli(datarun, cell, opts.triggers(trials), rasteropts, unmatched, 'hist_color', opts.hist_colors(p,:));

                if opts.plot_repr
                    hold on;
                    trials = repr(:,p) & singles;
                    [res{p,ncones+1}, ax2, histx, histy{p,ncones+1}] = rasterphli(datarun, cell, opts.triggers(trials), rasteropts, unmatched, 'hist_color', opts.hist_colors(p,:), 'MarkerStyle', 'o', 'LineStyle', '--');
                    if ~isempty(res{p,ncones+1}), numrast(p,ncones+1) = max(res{p,ncones+1}(:,2)); end
                    lax(p,ncones+1) = ax2(1);
                    rax(p,ncones+1) = ax2(2);

                    if opts.titles
                        title([num2str(p) ' stim/repr'])
                    end
                else
                    if opts.titles
                        title([num2str(p) ' stim'])
                    end
                end
            else
                trials = stim(:,p) & repr(:,q) & doubles;
                [res{p,q}, ax, histx, histy{p,q}] = rasterphli(datarun, cell, opts.triggers(trials), rasteropts, unmatched, 'hist_color', 'k');
                if opts.titles
                    title([num2str(p) ' stim, ' num2str(q) ' repr'])
                end
            end
            
            if ~isempty(res{p,q}), numrast(p,q) = max(res{p,q}(:,2)); end
            lax(p,q) = ax(1);
            rax(p,q) = ax(2);
        end 
    end
    
    if ~rasteropts.plot, continue; end
    
    if strcmp(opts.link_ylims, 'all')
        setaxesy(rax(:), 'Ylim', opts.set_ylims);
        set(rax(:,1:end-1), 'YTick', []);
    elseif strcmp(opts.link_ylims, 'rows')
        for p = 1:ncones
            setaxesy(rax(p,:), 'Ylim', opts.set_ylims);
            set(rax(p,1:end-1), 'YTick', []);
        end
    end
    
    for p = 1:ncones
        % Remove superfluous #trials ticks if all rasters have same number of trials
        lax_ylims = get(lax(p,:), 'YLim');
        if isequal(lax_ylims{:})
            set(lax(p,2:end), 'YTick', []);
        end
    end
    
    % Remove superfluous time axes
    set([rax(1:ncones-1,:); lax(1:ncones,:)], 'XTick', []);
end