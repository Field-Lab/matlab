function [rax,histx,histy,res] = additivity_raster_plot(datarun, stimstruct, cells, varargin)
% ADDITIVITY_RASTER_PLOT
% usage: [rax,histx,histy,res] = additivity_raster_plot(datarun, stimstruct, cells, varargin)
%
% 2011-08 phli
%

opts = inputParser;

% These are split out and passed to RASTER
opts.addParamValue('hist', true);
opts.addParamValue('hist_bin', 0.05);
opts.addParamValue('hist_bins', []);
opts.addParamValue('hist_line_width', 2);
opts.addParamValue('color', 0.5*[1 1 1]);
opts.addParamValue('MarkerSize', 5);
opts.addParamValue('triggers', datarun.triggers);
opts.addParamValue('start', -0.1);
opts.addParamValue('stop', []);
opts.addParamValue('plot', true);

% Unmatched will also be passed along to raster
opts.KeepUnmatched = true;

% These are for me
opts.addParamValue('hist_colors', {'r', [0 0.5 0], 'b'});
opts.addParamValue('intensities', []);
opts.addParamValue('link_ylims', false);
opts.addParamValue('with_blanks', true);
opts.addParamValue('sumhist_linestyle', 'k');
opts.addParamValue('subplot_width',  []);
opts.addParamValue('subplot_height', []);
opts.addParamValue('figure_by_cell_id', false);

opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;
if isempty(opts.stop)
    opts.stop = mean(diff(opts.triggers));
end;
rasteropts = keepfields(opts, 'hist', 'hist_bin', 'hist_bins', 'hist_line_width', 'color', 'MarkerSize', 'start', 'stop', 'plot');


if ~isempty(opts.intensities)
    intensities = opts.intensities;
else
    % Conventionally these stimuli have the same intensities for R, G, and B
    intensities = unique(stimstruct.rgbs(:,1));
    if ~opts.with_blanks
        intensities = setdiff(intensities, 0);
    end
    intensities = sort(intensities, 1, 'descend');
end

maps = stimstruct.mapindices;

numtrigg = length(opts.triggers);
if numtrigg < length(stimstruct.maps)
    warning(['Not enough triggers for all pulses.  Will only use the first ' num2str(numtrigg) ' pulses']);
end
if numtrigg > length(stimstruct.maps)
    error('More triggers than maps; I don''t know how to handle this.');
end

% Subplot dimensions
if isempty(opts.subplot_width),  opts.subplot_width  = length(maps);        end
if isempty(opts.subplot_height), opts.subplot_height = length(intensities); end

% New format is to leave cells as cellstruct that may have empties
if iscell(cells)
    cells = cell2mat(select(cells, @(e)(~isempty(e))));
end


for cell = cells
    if opts.plot
        if opts.figure_by_cell_id
            figure(cell);
        else
            figure;
        end
    end
    
    
    rax = [];
    for i = 1:length(intensities),
        intensity = intensities(i);
        lax = [];
        numrast = [];        
        
        % Plot row of rasters
        for j = 1:length(maps),
            map = maps(j);
            
            if opts.plot
                subplot(opts.subplot_height, opts.subplot_width, sub2ind([opts.subplot_width opts.subplot_height], j, i));
            end
            
            [res{i,j}, ax, histx, histy{i,j}] = rasterphli(datarun, cell, opts.triggers(stimstruct.rgbs(1:numtrigg,1) == intensity & (stimstruct.maps(1:numtrigg) == map)'), 'hist_color', opts.hist_colors{j}, rasteropts, unmatched);
            if ~isempty(res{i,j}), numrast(j) = max(res{i,j}(:,2)); end
            
            if ~opts.plot, continue; end
            
            lax(j)   = ax(1);
            rax(i,j) = ax(2);
            
            if (lax(j) > 0), set(lax(j), 'ytick', [], 'xtick', []); end
            
            if i == 1 && j == 1, title('1st'); end
            if i == 1 && j == 2, title('2nd'); end
            if i == 1 && j == 3, title('1+2'); end
            if j == 1, ylabel(num2str(intensity)); end
            if i ~= length(intensities) && rax(i,j), set(rax(i,j), 'xtick', 0); end
            if j ~= length(maps)        && rax(i,j), set(rax(i,j), 'yticklabel', []); end
        end
        if ~opts.plot, continue; end
        
        % Make stacked curve of hist from map0 and map1, plot to compare to
        % hist from map2 (the summed stimulus) (but only if there is a
        % summed stimulus response)
        if (length(lax) > 2 && lax(3) > 0)
            sumhist = histy{i,1} + histy{i,2};
            hold(rax(i,3), 'on');
            plot(rax(i,3), histx, sumhist, opts.sumhist_linestyle, 'LineWidth', opts.hist_line_width);
            maxsumhist = max(sumhist);
        else
            maxsumhist = 0;
        end

        
        % Link axes, saving right axis ylims first
        lax = lax(lax > 0);
        if ~isempty(lax), linkaxes(lax); end;
        for j = 1:size(rax, 2)
            ylims(j,:) = get(rax(i,j), 'ylim');
        end
        irax = rax(i,:);
        irax = irax(irax > 0);
        linkaxes(fliplr(irax));
        linkprop(fliplr(irax), 'YTickMode');
                
        % Set left axis ylims
        if ~isempty(lax) && ~isempty(numrast), set(lax(1), 'ylim', [0 max(numrast)]); end;
        
        % Adjust right axis ylims to include larger values if necessary, set to auto ticks
        ystart = min(ylims(:,1));
        if ystart < 0, ystart = 0; end
        yend   = max([ylims(:,2); maxsumhist]);
        set(rax(i,1), 'ylim', [ystart yend]);
        set(rax(i,1), 'YTickMode', 'auto');
        set(rax(i,1), 'xlim', [opts.start, opts.stop]);
    end
    
    % Normally we let each row have it's own y scaling (so that small
    % signals are visible).  Sometimes we want to link all y axes
    % though.
    if opts.link_ylims
        for i = 1:length(rax)
            ylims(i,:) = get(rax(i), 'ylim');
        end
        linkaxes(rax(:));
        set(rax(1), 'ylim', [min(ylims(:,1)) max(ylims(:,2))]);
    end

end
