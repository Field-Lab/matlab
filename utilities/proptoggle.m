function proptoggle(ax, prop)
% PROPTOGGLE    Toggle the property
% usage: newprop = proptoggle(ax, prop)
%
% 2011-07 phli

arrayfun(@(h) proptoggle1(h, prop), ax);

function proptoggle1(ax, prop)
set(ax, prop, toggle(get(ax, prop)));