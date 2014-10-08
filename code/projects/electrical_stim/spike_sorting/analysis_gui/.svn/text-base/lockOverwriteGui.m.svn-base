function response = lockOverwriteGui()

h.gui = figure('position', [500 500 200 100], 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'white');

uicontrol(h.gui,  'Style', 'text',...
    'String', 'Analysis is locked.  Old analysis will be overwritten.  Reanalyze anyway?',...
    'Position', [20 35 160 50], 'BackgroundColor', 'white');

h.yesButton = uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'yes', 'Position', [30 20 40 20], 'Callback', @yesFun);
h.noButton = uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'no', 'Position', [130 20 40 20], 'Callback', @noFun);


uiwait
%uiwait(h.gui)
close(h.gui)

    function yesFun(hObject, eventdata) %#ok<INUSD,INUSD>
        %close(h.gui)
        response = 1;
        uiresume
    end

    function noFun(hObject, eventdata) %#ok<INUSD,INUSD>
        %close(h.gui)
        response = 0;
        uiresume
    end
end