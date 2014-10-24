function plotTitles(dat,clr)

%-------------------------------------------
%
% plotTitles(dat)
%
% append titles to a plot
%
% dat must contain the fields:
% dat.cellType -- string with cell type (e.g. 'off midget')
% dat.dataSet -- string with data set (e.g. date/piece/run)
% dat.rgcId -- rgc identifier number
%
% optional input:
% clr -- color of text
%
% freeman, 3-19-2012
%-------------------------------------------


if ~exist('clr','var') || isempty(clr)
    clr = 'k';
end

% plot the cell type and id in the upper left corner
axes('Position',[0.025 0.975 0.2 0.2]);
axis off
h = text(0,0,(sprintf('%s, #%g',dat.cellType,dat.rgcId)));
set(h,'FontSize',20);
set(h,'Color',clr);

% plot the data set in the lower left corner
axes('Position',[0.025 0.025 0.2 0.2]);
axis off
h = text(0,0,lower(sprintf('%s',dat.dataSet)));
set(h,'FontSize',14);
set(h,'Color',clr);