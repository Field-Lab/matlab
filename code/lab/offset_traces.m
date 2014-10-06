function offset_traces(offset, ax)
% OFFSET_TRACES    Add incremental y-offsets to current plot traces
% usage: offset_traces(offset, ax)
%
% inputs: OFFSET    Amount in Y-Axis units to separate the traces by
%         AX        Axes handle, defaults to GCA
%
% This will try to keep track of the offset (in the 'offset' appdata for
% each trace).  So a call of OFFSET_TRACES(0) should restore the plot to
% normal.
%
% 2010-06 phli
%

if nargin < 2
    ax = gca;
end


children = get(ax, 'Children');
for i = 1:length(children)
    child = children(i);
    type = get(child, 'Type');
    
    if ~strcmp(type, 'line')
        continue;
    end
    
    this_offset = (i - 1) * offset;
    offset_trace(this_offset, child);
end




function offset_trace(offset, trace)
if ~(isempty(get(trace, 'XDataSource')) && isempty(get(trace, 'YDataSource')))
    warning('OFFSET_TRACES:LinkedData', 'Graph data may be linked to workspace variables; this will break offset tracking');
end

prev_offset = getappdata(trace, 'offset');
if isempty(prev_offset)
    prev_offset = 0;
end
delta = offset - prev_offset;

ydata = get(trace, 'YData');
ydata = ydata + delta;
set(trace, 'YData', ydata);

setappdata(trace, 'offset', offset);