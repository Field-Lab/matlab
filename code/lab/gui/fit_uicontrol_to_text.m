function fit_uicontrol_to_text(h, padding)
oldunits = get(h, 'Units');

if nargin < 2
    padding = guess_padding(h);
end

set(h, 'Units', 'characters');
ext = get(h, 'Extent');
pos = get(h, 'Position');
set(h, 'Position', [pos(1:2) ext(3:4)+padding]);

set(h, 'Units', oldunits);


% Doesn't seem like MathWorks is very consistent about how big the
% initially created uicontrol is, so we have to adhoc it.
function padding = guess_padding(h)
padding = [0 0];
style = get(h, 'Style');
if strcmp(style, 'pushbutton')
    padding = [3 1];
    return
end