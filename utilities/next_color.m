function nxt_color = next_color(ax, varargin)
% NEXT_COLOR    Decide what color to make the next plot onto held axes
% usage: nxt_color = next_color(ax, opts)
%
% inputs: ax         Axes handle; if not given or empty, defaults to gca
%
% opts:   color_order       []              The colors to be used.  NxRGB.
%
%         color_order_prop  'ColorOrder'    If color_order not set, then
%                                           looks in this axes property for colors
%
%         color_prop        'Color'         Looks at this property for each
%                                           plot element.  Can change to e.g.
%                                           MarkerFaceColor for markers.
%
%         Remaining opts are passed to findobj to limit which children of 
%         the axes are counted.  For example:
%            {'Type', 'line'}
%         If not given, then all children of the axes are used.
%
% Outputs the next color as a 1x3 RGB triple.
%
% Cycles through the existing lines on an axis and checks their colors.
% Counts how many times each color is used.  Picks the next color from the
% axes ColorOrder property that is least used.
%
% 2010-06 phli
%

if nargin < 1 || isempty(ax)
    ax = gca;
end

opts = inputParser;
opts.addParamValue('color_order',      []);
opts.addParamValue('color_order_prop', 'ColorOrder');
opts.addParamValue('color_prop',       'Color');
opts.KeepUnmatched = true;
opts.parse(varargin{:});
findopts = opts.Unmatched;
opts = opts.Results;



if isempty(fields(findopts))
    children = get(ax, 'Children');
else
    children = findobj('Parent', ax, findopts); % Luckily findobj is very tolerant and allows mixing of struct and paired vanilla opts.
end
num_children = length(children);



colors = opts.color_order;
if isempty(colors)
    colors = get(ax, opts.color_order_prop);
end
num_colors = size(colors, 1);


% Count up how many times each color has already been used within children
color_uses = zeros(num_colors, 1);
for i = 1:num_children
    child = children(i);
    child_color = get(child, opts.color_prop);

    for j = 1:num_colors
        color = colors(j, :);
        if all(color == child_color)
            color_uses(j) = color_uses(j) + 1;
        end
    end
end



least_num_uses = min(color_uses);
least_used_colors = find(color_uses == least_num_uses);
nxt_color = colors(least_used_colors(1), :);