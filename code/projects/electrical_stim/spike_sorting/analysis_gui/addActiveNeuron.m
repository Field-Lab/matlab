function [elecResp action] = addActiveNeuron(elecResp, movieIndex)


h.gui = figure('position', [500 500 200 150], 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'white');

uicontrol(h.gui, 'Style', 'text',...
    'String', 'Neuron to add to "active" list:',...
    'Position', [20 110 160 20], 'BackgroundColor', 'white');

h.neuronNo = uicontrol(h.gui, 'Style', 'edit', 'String', '', 'Position', [20 80 160 30],...
    'BackgroundColor', 'white', 'FontName', 'Courier');

h.yesButton = uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'add neuron', 'Position', [30 50 80 20], 'Callback', @addFun);
h.cancelButton = uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'cancel', 'Position', [130 50 40 20], 'Callback', @cancelFun);
h.createButton = uicontrol(h.gui,  'Style', 'pushbutton',...
    'String', 'create new ei', 'Position', [30 20 140 20], 'Callback', @createFun);

uiwait(h.gui)

close(h.gui)

    function addFun(hObject, eventdata) %#ok<INUSD,INUSD>
        newNeuron = str2double(get(h.neuronNo, 'String'));
        
        % check if specified neuron exists in params file
        datarun.names.rrs_params_path = elecResp.names.rrs_params_path;
        datarun = load_params(datarun, 'verbose', false, 'cell_type_depth', 2, 'sync_cell_ids', true);

        if ~any(datarun.cell_ids == newNeuron)
            warnH = warndlg('Neuron ID specified doesn''t exist in the params file');
            uiwait(warnH)
        else
            elecResp.cells.active{movieIndex} = [elecResp.cells.active{movieIndex} newNeuron];
            action = 'addExisting';
            
            elecResp.analysis.details.analysisFlags{movieIndex} = zeros(elecResp.stimInfo.nPulses(movieIndex),...
                length(elecResp.cells.active{movieIndex}));
            
            uiresume(h.gui)
        end
    end

    function cancelFun(hObject, eventdata) %#ok<INUSD,INUSD>
        action = 'cancel';
        uiresume(h.gui)
    end

    function createFun(hObject, eventdata)
        %set(h.gui, 'visible', 'off')
        ei = generateEIFromStimData(elecResp, movieIndex);
        if ~isempty(ei); %new ei successfully generated
            if max(elecResp.cells.all > 10000) %another generated ei already exists in this elecResp file
                newGenID = max(elecResp.cells.all) + 1;
            else
                newGenID = 10001;
            end
            elecResp.cells.allEIs{length(elecResp.cells.allEIs)+1} = ei;
            elecResp.cells.all = [elecResp.cells.all newGenID];
            elecResp.cells.active{movieIndex} = [elecResp.cells.active{movieIndex} newGenID];
            action = 'create';
        else
            action = 'cancel';
        end
        
        elecResp.analysis.details.analysisFlags{movieIndex} = zeros(elecResp.stimInfo.nPulses(movieIndex),...
            length(elecResp.cells.active{movieIndex}));
        
        uiresume(h.gui)

    end
end