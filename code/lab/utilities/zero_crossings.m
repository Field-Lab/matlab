function [output, polarity] = zero_crossings(traces, tracesx)
% ZERO_CROSSINGS    Find the zero crossing points (optionally polarity) of traces
% usage: [crossings, polarity] = zero_crossings(traces, tracesx)
%
% input: TRACES is a matrix of traces in columns
%        TRACESX (optional) is the x axis for the traces.
% 
% output: CROSSINGS is a cell array of the zero cross x points for each 
%         trace.  Each x-cross point is determined via linear interpolation
%         between the two points on either side of the crossing.
%
%         POLARITY if a cell array giving whether each crossing is 
%         positive-going (1) or negative going (-1)
%
% If no output is specified, the result will be plotted.
%
%
% 2010-06 phli
%

poses = traces >= 0;     % All positive points
crosses = diff(poses);   % 1 at positive crossings, -1 at negative crossings


% The rest we do trace by trace
num_traces  = size(traces,  2);
crossings = cell(num_traces, 1);
if nargout > 1
    polarity = cell(num_traces, 1);
end
for i = 1:num_traces
    trace = traces(:,i);
        
    if nargin < 2
        tracex = 1:length(trace);
    elseif size(tracesx, 1) == 1
        tracex = tracesx;
    else
        tracex = tracesx(i,:);
    end
    
    % Calculate the zero crossing x point via linear interpolation
        precrosses = find(crosses(:,i) ~= 0);
        postcrosses = precrosses + 1;
        
        y1 = trace(precrosses);
        y2 = trace(postcrosses);
        deltaxzerofactor = y1 ./ (y1 - y2);
        
        x1 = tracex(precrosses);
        x2 = tracex(postcrosses);
        deltax = x2 - x1;
        crossings{i} = x1 + (deltax .* deltaxzerofactor');
    % --
    
    if nargout > 1
        polarity{i} = crosses(crosses(:,i) ~= 0, i)';
    end    
end


% Plot, if no other output
if nargout < 1 && num_traces == 1
    plot(tracex, trace);
    hold on;
    plot(crossings{1}, 0, 'ok');
end

if nargout > 0
    output = crossings;
end