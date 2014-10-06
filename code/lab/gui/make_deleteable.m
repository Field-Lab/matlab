function make_deleteable(h)

% Get element's context menu
cmenu = get(h, 'UIContextMenu');

% Ensure element has a context menu
if isempty(cmenu)
    cmenu = uicontextmenu;
    set(h, 'UIContextMenu', cmenu);
end

% Add delete to context menu
uimenu(cmenu, 'Label', 'Delete', 'Callback', @(varargin)(delete(h)));