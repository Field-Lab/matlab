function [eipicker csdpicker] = add_ei_csd_toggle(fig, datarun, cellid, position)

% FIXME: Put these in the app in case PICKCSD below needs to call GET_CSD.
% GET_CSD shouldn't really depend on datarun; it should be able to run with
% just the EI and the POSITIONS, but getting the neighbor_struct is 
% currently pretty kludgy so do it this way for now.
setappdata(fig, 'datarun', datarun);
setappdata(fig, 'cellid', cellid);

% Buttons to choose between EI/CSD display
eipicker = uicontrol(fig, ...
    'Style'   , 'pushbutton',      ...
    'String'  , 'Switch to CSD',   ...
    'FontSize', 12,                ...
    'Units'   , 'normalized',      ...
    'Callback', @pickcsd           ...
);
if nargin > 3, 
    if length(position) == 2, position(3:4) = [1 1]; end
    set(eipicker, 'Position', position); 
end
fit_uicontrol_to_text(eipicker);
setappdata(fig, 'eipicker', eipicker);

csdpicker = uicontrol(fig, ...
    'Style'   , 'pushbutton',              ...
    'Visible' , 'off',                     ...
    'String'  , 'Switch to EI',            ...
    'FontSize', get(eipicker, 'FontSize'), ...
    'Units'   , get(eipicker, 'Units'),    ...
    'Position', get(eipicker, 'Position'), ...
    'Callback', @pickei                    ...
);
setappdata(fig, 'csdpicker', csdpicker);


function pickcsd(h,~)
f = getfig(h);

if isempty(getappdata(f, 'csd')),
    oldstr = get(h, 'String');
    set(h, 'String', 'Loading...');
    drawnow
    
    % FIXME: See ADD_EI_CSD_TOGGLE above
    datarun = getappdata(f, 'datarun');
    cellid = getappdata(f, 'cellid');
    csd = get_csd(datarun, cellid);
    setappdata(f, 'csd', csd);
    
    set(h, 'String', oldstr);
end
setappdata(f, 'data', getappdata(f, 'csd'));

% Switch buttons
set(getappdata(f, 'csdpicker'), 'Visible', 'on');
set(h, 'Visible', 'off');

replot(f);



function pickei(h,~)
f = getfig(h);

setappdata(f, 'data', getappdata(f, 'ei'));

% Switch buttons
set(getappdata(f, 'eipicker'), 'Visible', 'on');
set(h, 'Visible', 'off');

% Replot
replot(f);


function replot(f)
api = getappdata(f, 'api');
slider = getappdata(f, 'ei_slider');
api.slider_plot(slider, []);