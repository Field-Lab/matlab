function clusterBin = manualClustering(traces, channels)

%Traces - 3-dimensional array: traces*electrodes*samples.
%Channels - vector of electrode numbers corresponding to 2nd dimension of traces

%% Construct the gui components

% initalizing variables
Data.placeHolder = 1;
i = 1;
j = 1;

Data.electrodes = channels;
nDisplayElectrodes = length(channels);

Data.traces = traces;

H.mainGui = figure('Visible', 'off', 'Position', [50 100 1200 600], 'Name', 'Hand-Clustering Main.', 'Resize', 'off');

%generates axes and associated pushbuttons
axesSpacing = (1100 - (nDisplayElectrodes-1)*5)/nDisplayElectrodes;
H.axes = cell(1, nDisplayElectrodes);
H.pushbuttons = cell(1, nDisplayElectrodes);

for i = 1:nDisplayElectrodes %#ok<FXUP>
    H.axes{i} = axes('Parent', H.mainGui, 'Units', 'pixels', 'Position', [50 + round((axesSpacing+5)*(i-1)) 50 round(axesSpacing) 200]);
    H.pushbuttons{i} = uicontrol(H.mainGui, 'Style', 'pushbutton', 'String', ['electrode ' num2str(Data.electrodes(i))],...
        'Position', [50 + round((axesSpacing+5)*(i-1)) 250 round(axesSpacing) 30], 'Callback', {@electrodeCB, i});
end
guidata(H.mainGui, H)


%displays data

% if statements are for the situation where a button is pressed before all figures have loaded
traceColor = cool(size(Data.traces,1));
for i = 1:nDisplayElectrodes           %#ok<FXUP>
    if sum(allchild(0) == H.mainGui)
        axes(H.axes{i})
    end
    if sum(allchild(0) == H.mainGui)
        set(H.axes{i},'YTickLabel',{''})
        hold on
        for j = 1:size(Data.traces,1) %#ok<FXUP>
            current = plot(squeeze(Data.traces(j,i,:)));
            set(findobj(current,'Type','line'), 'Color', traceColor(j,:))
        end
        hold off
    end
end



%% initialization tasks

if sum(allchild(0) == H.mainGui)
    set(H.mainGui, 'Visible', 'on','Color', [0.7 0.8 1])
    guidata(H.mainGui, H)
end

%% Callbacks for gui

    function electrodeCB(hObject, eventdata, index) %#ok<INUSL>
        H = guidata(H.mainGui);
        H.electrodeChoice = index;
        guidata(H.mainGui, H); %probably unecessary
        Data.clusterBin = manualCluster2D(squeeze(Data.traces(:,H.electrodeChoice,:)));
        guidata(H.mainGui,Data); %so that clusterBin can be returned
        close(H.mainGui);
    end


%% other

if sum(allchild(0) == H.mainGui)
    uiwait(H.mainGui)
end

%sets return value when gui figure is closed by function electrodeCB
clusterBin = Data.clusterBin;

end %end of main function
