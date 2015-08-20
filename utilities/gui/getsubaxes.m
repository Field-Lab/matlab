function ax = getsubaxes(fig)

if nargin < 1, fig = gcf; end
children = get(fig, 'Children');
ax = findobj(children, 'flat', 'Type', 'axes');