function outstructs = raster_subplot(instructs, varargin)
% RASTER_SUBPLOT    Flexibly plot rasters into a subplotted panel
% usage: outstructs = raster_subplot(instructs, opts)
% 
% INSTRUCTS is a cell array of struct arrays.  The dimensions of the cell
% array are the layout of subplots.  Each struct array can have multiple
% entries if multiple PSTHs should be plotted on the same subplot.
%
% Each struct should have fields:
%   datarun
%   cellspec
%   triggers
%   
% Each struct optionally can include raster opts:
%   hist_bin
%   start
%   stop
%   color
%   hist_color
%   hist_line_width
%   LineStyle
%   MarkerStyle
%   MarkerSize
%
% 2012 phli
%

opts = inputParser();
opts.addParamValue('cla', true);
opts.addParamValue('add_width' , 0);
opts.addParamValue('add_height', 0);
opts.addParamValue('raster_offset', [0 0]);
opts.addParamValue('setaxesy', true);
opts.addParamValue('YLim', true);
opts.addParamValue('titles', {});
opts.addParamValue('shadex', []);
opts.addParamValue('pbaspect', []);
opts.parse(varargin{:});
opts = opts.Results;

rasters_size = size(instructs);
subplot_height = rasters_size(1) + opts.add_height;
subplot_width  = rasters_size(2) + opts.add_width;

for i = 1:rasters_size(1)
    for j = 1:rasters_size(2)
        ins = instructs{i,j};
        if isempty(ins), continue; end
        
        sanesubplot(subplot_height, subplot_width, {i+opts.raster_offset(1) j+opts.raster_offset(2)});
        outstructs{i,j} = raster_subplot_(ins, opts);
        
        if hasindex(opts.titles, [i j]), title(opts.titles{i,j}); end
    end
end


ax = cell2mat(cellfun(@(s)(vertcat(s.AX)), invertselect(outstructs(:), @isempty), 'UniformOutput', false));
lax = ax(:,1);
rax = ax(:,2);
rax = rax(rax > 0 & lax > 0);
lax = lax(lax > 0);

if opts.setaxesy, setaxesy(rax, 'YLim', opts.YLim); end



function outstructs = raster_subplot_(instructs, opts)
hold on;
opts.rastax = gca;
opts.histax = 0;

for i = 1:length(instructs)
    outstructs(i) = raster_subplot__(instructs(i), opts);
    opts.histax = outstructs(i).AX(2);

    if i == 1,
        axes(opts.histax);
        hold on;
    end
end

axes(opts.rastax);
hold off
axes(opts.histax);
hold off


function outstruct = raster_subplot__(instruct, opts)
instruct.opts.plot = true;
instruct.opts.hist = true;
instruct.opts.rastax = opts.rastax;
instruct.opts.histax = opts.histax;
instruct.opts.shadex = opts.shadex;
instruct.opts.pbaspect = opts.pbaspect;

[outstruct.res,outstruct.AX,outstruct.histx,outstruct.tt] = rasterphli(instruct.datarun, instruct.cellspec, instruct.triggers, instruct.opts);