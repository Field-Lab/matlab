function [newNeuronID action] = chooseNeuronID(params_path)

h.chooseIDGui = figure('position', [500 500 200 150], 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'white');

uicontrol(h.chooseIDGui, 'Style', 'text',...
    'String', 'Choose an ID number for the new EI:',...
    'Position', [20 90 160 40], 'BackgroundColor', 'white');

h.neuronNo = uicontrol(h.chooseIDGui, 'Style', 'edit', 'String', '', 'Position', [20 60 160 30],...
    'BackgroundColor', 'white', 'FontName', 'Courier');

h.yesButton = uicontrol(h.chooseIDGui,  'Style', 'pushbutton',...
    'String', 'add neuron', 'Position', [30 30 80 20], 'Callback', @addFun);
h.cancelButton = uicontrol(h.chooseIDGui,  'Style', 'pushbutton',...
    'String', 'cancel', 'Position', [130 30 40 20], 'Callback', @cancelFun);

uiwait(h.chooseIDGui)

close(h.chooseIDGui)


    function addFun(hObject, eventdata) %#ok<INUSD,INUSD>
        newNeuronID = str2double(get(h.neuronNo, 'String'));

        % check if specified neuron exists in params file
        datarun.names.rrs_params_path = params_path;
        datarun = load_params(datarun, 'verbose', false, 'cell_type_depth', 2, 'sync_cell_ids', true);

        if any(datarun.cell_ids == newNeuronID)
            warnH = warndlg(['Neuron ID specified already exists in ' params_path]);
            uiwait(warnH)
        else
            action = 'success';
            uiresume(h.chooseIDGui)
        end
    end

    function cancelFun(hObject, eventdata) %#ok<INUSD,INUSD>
        newNeuronID = 0;
        action = 'cancel';
        uiresume(h.chooseIDGui)
    end
end
