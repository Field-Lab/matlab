function isAligned = checkAlignment(traces, traceMins)

h.gui = figure('position', [500 400 600 400], 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'white',...
    'visible', 'off', 'name', 'check alignment');

uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'confirm alignment', 'Position', [50 20 200 20], 'Callback', @confirm);
uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'cancel', 'Position', [350 20 200 20], 'Callback', @cancel);

axes('parent', h.gui, 'units', 'pixels', 'position', [50 60 500 320]);

hold on
for i = 1:length(traceMins)
    plot(traces(i,traceMins(i)-10:traceMins(i)+15))
end
hold off

set(h.gui, 'visible', 'on')

uiwait(h.gui)

    function confirm(hObject, eventdata) %#ok<INUSD>
        isAligned = 1;
        close(h.gui)
    end

    function cancel(hObject, eventdata) %#ok<INUSD>
        isAligned = 0;
        close(h.gui)
    end
end