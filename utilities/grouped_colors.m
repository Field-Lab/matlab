function colors = grouped_colors(groups, num, outgroup_color, colors_or_axes)
% GROUPED_COLORS    Create a matrix of RGB triples for the groups given
% usage: colors = grouped_colors(groups, [num, outgroup_color, colors_or_axes])
%
% inputs: groups   Cell array of M group vectors (see below)
%
%         num      The total number of RGB triples to output.  Defaults to
%                  the highest index given in any of the group vectors
%
%         outgroup_color   The color to give to all indices not included in
%                          a group.  Defaults to the next available color
%                          in the colorset (see below) or else to white.
%
%         colors_or_axes   Either an Nx3 matrix of RGB triples to use, N >=
%                          M, or else an axes handle to get the ColorOrder
%                          from.  Defaults to gca.
%
% The groups cell array should contain N vectors, each made up of a series
% of indices that should be given the same color.  E.g.: {[1 3 5] [2 4 6]},
% where this is coloring the odd indices one color and the even indices a
% second color.  If num was set to 10, there would also be a third color
% assigned to indices 7:10.
%
% 2010-06 phli
%

if nargin < 4
    colors_or_axes = gca;
end


if isscalar(colors_or_axes)
    % It's a axes handle
    colorset = get(colors_or_axes, 'ColorOrder');
else
    % It's RGB triples
    colorset = colors_or_axes;
end


num_groups = length(groups);
if nargin < 3
    if size(colorset, 1) > num_groups
        outgroup_color = colorset(num_groups+1, :);
    else
        outgroup_color = [1 1 1]; % Make it white I guess...
    end
end


if nargin < 2
	maxes = collect(groups(:), @(group) max(group(:)));
    num = max(cell2mat(maxes));
end


colors = repmat(outgroup_color, num, 1);
for i = 1:length(groups)
    group = groups{i};
    color = colorset(i, :);
    colors(group(:),:) = repmat(color, length(group), 1);
end