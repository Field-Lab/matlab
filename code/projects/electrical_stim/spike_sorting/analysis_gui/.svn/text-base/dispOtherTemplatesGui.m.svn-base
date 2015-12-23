function [ otherDispEIs ] = dispOtherTemplatesGui(main)


chooseTemplatesToDispGui = figure('position', [500 500 300 250], 'color', 'white', 'Toolbar', 'none',...
    'Menubar', 'none');

uicontrol(chooseTemplatesToDispGui, 'Style', 'text', 'String', 'neuron id: ',...
    'Position', [20 220 180 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
uicontrol(chooseTemplatesToDispGui, 'Style', 'text', 'String', 'external EI? ',...
    'Position', [200 220 90 15], 'HorizontalAlignment', 'right', 'BackgroundColor', 'white');

neuronText{1} = uicontrol(chooseTemplatesToDispGui, 'Style', 'edit', 'String', '', 'Position', [20 180 180 25]);
neuronText{2} = uicontrol(chooseTemplatesToDispGui, 'Style', 'edit', 'String', '', 'Position', [20 150 180 25]);
neuronText{3} = uicontrol(chooseTemplatesToDispGui, 'Style', 'edit', 'String', '', 'Position', [20 120 180 25]);
neuronText{4} = uicontrol(chooseTemplatesToDispGui, 'Style', 'edit', 'String', '', 'Position', [20 90  180 25]);
neuronText{5} = uicontrol(chooseTemplatesToDispGui, 'Style', 'edit', 'String', '', 'Position', [20 60  180 25]);


externalCB{1} = uicontrol(chooseTemplatesToDispGui, 'Style', 'checkbox', 'position', [260 185 20 15], 'string', '', 'Value', 0);
externalCB{2} = uicontrol(chooseTemplatesToDispGui, 'Style', 'checkbox', 'position', [260 155 20 15], 'string', '', 'Value', 0);
externalCB{3} = uicontrol(chooseTemplatesToDispGui, 'Style', 'checkbox', 'position', [260 125 20 15], 'string', '', 'Value', 0);
externalCB{4} = uicontrol(chooseTemplatesToDispGui, 'Style', 'checkbox', 'position', [260 95  20 15], 'string', '', 'Value', 0);
externalCB{5} = uicontrol(chooseTemplatesToDispGui, 'Style', 'checkbox', 'position', [260 65  20 15], 'string', '', 'Value', 0);

uicontrol(chooseTemplatesToDispGui,  'Style', 'pushbutton', 'String', 'accept',...
    'Position', [20 20 120 20], 'Callback', @acceptFunction);
uicontrol(chooseTemplatesToDispGui,  'Style', 'pushbutton', 'String', 'cancel',...
    'Position', [160 20 120 20], 'Callback', @cancelFunction);

uiwait(chooseTemplatesToDispGui)


    function acceptFunction(~, ~)
        
        count = 0;
        for ii = 1:5
            neuronString = get(neuronText{ii}, 'string');
            if ~isempty(neuronString)
                count = count + 1;
                filePath = get(main.filePath, 'String');
                if get(externalCB{ii}, 'value')
                    if exist([filePath filesep 'eiFile' neuronString '.mat'], 'file')
                        eiData = load([filePath filesep 'eiFile' neuronString '.mat']);
                        otherDispEIs{count} = eiData.ei;
                    else
                        error(['No external EI file found for neuron ' neuronString])
                    end
                else
                    fileName = get(main.fileName, 'String');
                    temp = load([filePath fileName]);
                    elecResp = temp.elecResp;
                    cellInd = find(elecResp.cells.all==str2double(neuronString));
                    if ~isempty(cellInd)
                        otherDispEIs{count} = elecResp.cells.allEIs{cellInd};
                    else
                        error(['No EI found for neuron ' neuronString...
                            ' in elecResp file -- neuron doesn''t exist or doesn''t have large enough signal (>=10 DAQ) on the recording electrode(s).'])
                    end
                end
            end
        end
        if count == 0
            otherDispEIs = [];
        end
        close(chooseTemplatesToDispGui)
    end

    function cancelFunction(~, ~)
        otherDispEIs = 0; %tells main GUI to not change current values for main.otherDispEIs
        close(chooseTemplatesToDispGui)
    end
end

