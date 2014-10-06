function findNeuronsOnElecsGui(elecResp, main)



%% Construct the gui components

%electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%channelsToUse = electrodeMap.getAdjacentsTo(elecResp.cells.recElec, 1)';
if isfield(elecResp.cells,'mainEI')
if size(elecResp.cells.mainEI,1) == 512 %LG modification 1/24/2014
    channelsToUse = getCluster512(elecResp.cells.recElec);
else
    channelsToUse = getCluster(elecResp.cells.recElec);
end
else
    channelsToUse = getCluster512(elecResp.cells.recElec);
    disp('Assumes 512-electrode array'); 
%     channelsToUse = getCluster(elecResp.cells.recElec);
%     disp('Assumes 61-electrode array');
end

nChannels = length(channelsToUse);

mainPos = get(main.gui, 'Position');
panelPos = get(main.bottomPanel, 'position');
axesPos = cell(1, nChannels);
for i = 1:nChannels
    axesPos{i} = get(main.aEi{1,i}, 'position');
end


guiPos = [mainPos(1)+panelPos(1)+axesPos{1}(1) - 25, mainPos(2)+panelPos(2)+axesPos{1}(2) + axesPos{1}(4)+10,...
    axesPos{nChannels}(1) + axesPos{nChannels}(3) - axesPos{1}(1) + 40, axesPos{1}(4)+200];
offsetFromPanel = axesPos{1}(1) - 25;



%% construct the gui

findNeuronsGui = figure('position', [500 500 300 180], 'color', 'white', 'Toolbar', 'none',...
    'Menubar', 'none');

uicontrol(findNeuronsGui, 'Style', 'text', 'String', 'electrodes to look for signal on: ',...
    'Position', [20 140 180 15], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

channelCheckBox = cell(length(channelsToUse));
for i = 1:length(channelsToUse)
    if i <= 6
        channelCheckBox{i} = uicontrol(findNeuronsGui, 'Style', 'checkbox', 'Units', 'pixels',...
            'Position', [20+45*(i-1) 110 40 20], 'String', num2str(channelsToUse(i)),...
            'Fontweight', 'bold', 'Value', 0);
    elseif i <= 12
        channelCheckBox{i} = uicontrol(findNeuronsGui, 'Style', 'checkbox', 'Units', 'pixels',...
            'Position', [20+45*(i-7) 90 40 20], 'String', num2str(channelsToUse(i)),...
            'Fontweight', 'bold', 'Value', 0);
    else
        channelCheckBox{i} = uicontrol(findNeuronsGui, 'Style', 'checkbox', 'Units', 'pixels',...
            'Position', [20+45*(i-13) 70 40 20], 'String', num2str(channelsToUse(i)),...
            'Fontweight', 'bold', 'Value', 0);
    end
end

uicontrol(findNeuronsGui, 'Style', 'text', 'String', 'threshold (DAQ units; min = 10): ',...
    'Position', [20 50 180 15], 'HorizontalAlignment', 'right', 'BackgroundColor', 'white');
thresholdText = uicontrol(findNeuronsGui, 'Style', 'edit', 'String', '20',...
    'Position', [210 50 30 25], 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');

uicontrol(findNeuronsGui,  'Style', 'pushbutton', 'String', 'find neurons',...
    'Position', [20 20 120 20], 'Callback', @findNeuronsDummy);
uicontrol(findNeuronsGui,  'Style', 'pushbutton', 'String', 'cancel',...
    'Position', [160 20 120 20], 'Callback', @cancelFunction);

uiwait(findNeuronsGui)


    function findNeuronsDummy(hObject, eventdata) %#ok<INUSD,INUSD>
        thresh = str2double(get(thresholdText, 'String'));
        channelsOriginal = channelsToUse;
        for j = length(channelsToUse):-1:1
            if ~get(channelCheckBox{j}, 'Value')
                channelsToUse(j) = [];
            end
        end
        nChannels = length(channelsOriginal);
        if isfield(elecResp.cells, 'all')
            neurons = findNeuronsAboveThresh(elecResp.names.rrs_ei_path, elecResp.names.rrs_params_path,...
                channelsToUse, thresh, elecResp.cells.all);
        else
            neurons = findNeuronsAboveThresh(elecResp.names.rrs_ei_path, elecResp.names.rrs_params_path,...
                channelsToUse, thresh);
        end
        close(findNeuronsGui)
        
        % make figure with eis of neurons with signal above threshold
        if ~isempty(neurons)
            eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(elecResp.names.rrs_ei_path);
            eiData = cell(length(neurons), 1); %stores EIs of neurons: eiData{neuronIndex}(channels, samples)
            eiYMin = 5000; eiYMax = -5000;
            eiMinPos = zeros(length(neurons), 1);
            for j = 1:length(neurons)
                eiData{j} = eiFile.getImage(neurons(j));
                eiData{j} = squeeze(eiData{j}(1, 2:end, :));
                eiYMin = min([min(min(eiData{j})), eiYMin]);
                eiYMax = max([max(max(eiData{j})), eiYMax]);
                % if main recording electrode is one of the searched electrodes, use it to determine
                % template minimum
                if any(channelsToUse == elecResp.cells.recElec)
                    eiMinPos(j) = find(squeeze(eiData{j}(elecResp.cells.recElec, :))...
                        ==min(squeeze(eiData{j}(elecResp.cells.recElec,:))), 1);
                else %otherwise use first electrode in list of chosen to determine template minimum
                    eiMinPos(j) = find(squeeze(eiData{j}(channelsToUse(1), :))...
                        ==min(squeeze(eiData{j}(channelsToUse(1), :))), 1);
                end
            end
            clear eiFile
            

            neuronsOnElecsFig = figure('Position', guiPos, 'color', 'white', 'Menubar', 'none');

            neuronColors = hsv(length(neurons));
            for j = 1:nChannels


                %axes('Parent', neuronsOnElecsFig, 'Units', 'pixels', 'Position',...
                %    [axesPos{j}(1) - offsetFromPanel, 20, axesPos{j}(3), axesPos{j}(4)])
                
                axes('Parent', neuronsOnElecsFig, 'Position',...
                    [(axesPos{j}(1) - offsetFromPanel)/guiPos(3), 20/guiPos(4), axesPos{j}(3)/guiPos(3), axesPos{j}(4)/guiPos(4)])

 
                hold on
                for k = 1:length(neurons)
                    startInd = eiMinPos(k)-10;
                    if startInd < 1
                        plot(2-startInd:eiMinPos(k)+16-startInd, eiData{k}(channelsOriginal(j), 1: eiMinPos(k) + 15), 'LineWidth', 2, 'Color', neuronColors(k,:));
%                         dummyTrace(j,k,:) = eiData{k}(channelsOriginal(j), 1: eiMinPos(k) + 15); 
                    else
                        plot(eiData{k}(channelsOriginal(j), eiMinPos(k) - 10: eiMinPos(k) + 15), 'LineWidth', 2, 'Color', neuronColors(k,:));
%                         dummyTrace(j,k,:) = eiData{k}(channelsOriginal(j), eiMinPos(k) - 10: eiMinPos(k) + 15);
                    end
                end
                 
                set(gca, 'YLim', [eiYMin eiYMax]);
                title(['electrode ' num2str(channelsOriginal(j))], 'FontSize', 14, 'FontWeight', 'bold')

            end

            %legend
            axes('parent', neuronsOnElecsFig,...
                'position', [(axesPos{1}(1) - offsetFromPanel)/guiPos(3) (axesPos{j}(4)+50)/guiPos(4) 200/guiPos(3) 140/guiPos(4)],...
                'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1])
            hold on
            for k = 1:length(neurons)
                plot(0, 0.9*k/length(neurons)+0.05, 's', 'MarkerEdgeColor', neuronColors(k,:),...
                    'MarkerFaceColor', neuronColors(k,:), 'MarkerSize', 10)
                text(0.1, 0.9*k/length(neurons)+0.05, ['neuron ' num2str(neurons(k))])
            end
            hold off
            
        else
            warndlg('No neurons found with signal above specified threshold on specified electrode(s).')
            
        end
    end

    function cancelFunction(hObject, eventdata) %#ok<INUSD,INUSD>
        close(findNeuronsGui)
    end


end