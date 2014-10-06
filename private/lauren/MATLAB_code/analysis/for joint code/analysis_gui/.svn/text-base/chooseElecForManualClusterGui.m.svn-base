function chosenElec = chooseElecForManualClusterGui(main, channels)

%% preparing for gui construction

nChannels = length(channels);

mainPos = get(main.gui, 'Position');
panelPos = get(main.bottomPanel, 'position');

axesPos = cell(1, nChannels);
for i = 1:nChannels
    axesPos{i} = get(main.aEi{1,i}, 'position');
end

guiPos = [mainPos(1)+panelPos(1)+axesPos{1}(1) - 10, mainPos(2)+panelPos(2)+axesPos{1}(2) - 120,...
    axesPos{nChannels}(1) + axesPos{nChannels}(3) - axesPos{1}(1) + 20, 60];

offsetFromPanel = axesPos{1}(1) - 10;

%% Construct the gui components

chooseElec.gui = figure('Visible', 'off', 'Position', guiPos, 'color', 'white',...
    'Toolbar', 'none', 'Menubar', 'none');

%generates axes and associated pushbuttons
chooseElec.pushButtons = cell(1, nChannels);

uicontrol(chooseElec.gui, 'Style', 'text', 'String', 'choose an electrode for manual clustering',...
    'Position', [round(guiPos(3)/2)-150 40 300 15], 'fontsize', 12,...
    'fontweight', 'bold', 'BackgroundColor', 'white');


for i = 1:nChannels
    chooseElec.pushButtons{i} = uicontrol(chooseElec.gui, 'Style', 'pushbutton',...
        'String', ['electrode ' num2str(channels(i))],...
        'Position', [axesPos{i}(1) - offsetFromPanel, 10, axesPos{i}(3), 20],...
        'Callback', {@electrodeCB, i});
end


%% initialization tasks

%if sum(allchild(0) == H.mainGui)
set(chooseElec.gui, 'Visible', 'on')
%end

uiwait(chooseElec.gui)


%% Callbacks for gui

    function electrodeCB(hObject, eventdata, index) %#ok<INUSL>
        chosenElec = channels(index);
        close(chooseElec.gui);
    end
end