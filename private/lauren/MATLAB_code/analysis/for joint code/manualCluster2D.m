function clusterBin = manualCluster2D(traces)
% manualCluster2D: manual clustering of waveforms using a gui
% 
% manualCluster2D(traces) displays two groups of traces, then allows the user to move traces from 
% group to group by drawing a polygon around particular traces themselves or their representation in
% PC space.
%
% arguments
%   traces: set of traces on a single electrode or multiple electrodes
%   concatenated together(traces x samples)
%
% returns
%   clusterBin: a vector (length = number of traces) of ones and twos, denoting which traces
%   are in the category "containing target spikes" (left plot, denoted by a one) or "other traces"
%   (right plot, denoted by a zero)
%
% requirements:
% 	icons.mat to display correctly
% 	moveTraces, lassoTraces.m (used by moveTraces), lasso.m (used by moveTraces), PCChooser.m (used
%   by moveTraces), and updateManualClusterAxes to function correctly

% written by Lauren Hruby, SNL-E
% 2008-10-24
% edited 2008-11-25: changed clusterBin to be 1 for artifact and 2 for spikes (vs. 0 and 1)
% bug: first attempt at moving traces fails (noticed 2008-11-25)
% bug fixed 2008-12-16
% branched 2009-07-02 for use with new analysis Gui
%   changes: changed clusterIndeces to clusterBin, and changed 1 vs. 2 to 0 vs. 1 in clusterBin

Data.traces = traces;

Icons.initializer = 1;
try
    load('icons.mat');
catch
end
    
i = 1; %#ok<NASGU> % so that i can be used later in static workspace

Data.nTraces  = size(Data.traces, 1);
Data.nSamples = size(Data.traces, 2);

% calcutes values to ues for axes limits, rounding to 50 after adding a buffer of 20 on each side
Data.yMin = 50*floor((min(min(Data.traces)) - 20)/50);
Data.yMax = 50*ceil ((max(max(Data.traces)) + 20)/50);

Data.includedIndeces = [];
Data.discludedIndeces = 1:Data.nTraces;
Data.includedTraces  = [];
Data.discludedTraces = Data.traces;

Data.prevDiscludedIndeces = Data.discludedIndeces;
Data.prevIncludedIndeces  = Data.includedIndeces;

H.singElecGui = figure('Visible', 'off', 'Position', [50 100 1000 700], 'Name', 'Hand-Clustering Single-Electrode.', 'Resize', 'off');
guidata(H.singElecGui, [])% clears any old guidata



%% Construct the gui components

% axes
H.axesLeft = axes('Parent', H.singElecGui, 'Units', 'pixels', 'Position', [50 50 400 400]);
H.axesRight = axes('Parent', H.singElecGui, 'Units', 'pixels', 'Position', [550 50 400 400]);

% buttons
H.selectionType =   uibuttongroup('Parent', H.singElecGui, 'Visible', 'off', 'Units', 'pixels',...
	 'BorderType', 'none',...
     'BackgroundColor', [0.7 0.8 1],'Position', [50 600 200 60]);
H.usePCAButton =    uicontrol('Style', 'Radio', 'String', 'select using PCA plot',...
     'BackgroundColor', [0.7 0.8 1], 'Pos', [0 30 200 30], 'Parent', H.selectionType);
H.useTracesButton = uicontrol('Style', 'Radio', 'String', 'select using traces',...
     'BackgroundColor', [0.7 0.8 1], 'Pos', [0 0 200 30], 'Parent', H.selectionType);
 
H.undoButton = uicontrol(H.singElecGui,   'Style', 'pushbutton', 'String', 'undo',...
    'Position', [550 630 100 15]);
H.acceptButton = uicontrol(H.singElecGui, 'Style', 'pushbutton', 'String', 'accept clusters',...
    'Position', [550 600 100 15]);

try
    H.toLeftButton = uicontrol(H.singElecGui,  'Style', 'pushbutton', 'CData', Icons.leftArrow,...
        'Position', [480 300 24 24]);
    H.toRightButton = uicontrol(H.singElecGui, 'Style', 'pushbutton', 'CData', Icons.rightArrow,...
        'Position', [480 250 24 24]);
    H.swapButton    = uicontrol(H.singElecGui, 'Style', 'pushbutton', 'CData', Icons.swapArrow,...
        'Position', [480 200 24 24]);
catch % if icons.mat wasn't loaded correctly
    H.toLeftButton = uicontrol(H.singElecGui,  'Style', 'pushbutton', 'String', '<--',...
        'Position', [480 300 24 24]);
    H.toRightButton = uicontrol(H.singElecGui, 'Style', 'pushbutton', 'String', '-->',...
        'Position', [480 250 24 24]);
    H.swapButton    = uicontrol(H.singElecGui, 'Style', 'pushbutton', 'String', '<-->',...
        'Position', [480 200 24 24]);
end
 


% so that H is fully formed before used as an argument:
set(H.undoButton,   'Callback', {@undoPreviousMove})
set(H.acceptButton, 'Callback', {@acceptClusters})

set(H.toLeftButton,  'Callback', {@moveTracesDummy, H})
set(H.toRightButton, 'Callback', {@moveTracesDummy, H})
set(H.swapButton,    'Callback', {@moveTracesDummy, H})

% labels for graphs
uicontrol(H.singElecGui, 'Style', 'text', 'String', 'traces containing target spikes',...
    'Fontsize', 14,...
    'HorizontalAlignment', 'left', 'BackgroundColor', [0.7 0.8 1], 'Position', [50 460 400 30])
uicontrol(H.singElecGui, 'Style', 'text', 'String', 'other traces',...
    'Fontsize', 14,...
    'HorizontalAlignment', 'left', 'BackgroundColor', [0.7 0.8 1], 'Position', [550 460 400 30])

%% initialization tasks

set(H.singElecGui, 'Visible', 'on', 'Color', [0.7 0.8 1])

% Initializes button group "selectionType" properties. 
set(H.selectionType, 'SelectedObject', [], 'Visible', 'on');

% displays data
set(H.axesLeft,  'Ylim', [Data.yMin Data.yMax])
set(H.axesRight, 'YLim', [Data.yMin Data.yMax])

axes(H.axesRight)
H.traceColor = cool(Data.nTraces);
hold on
for i = 1:Data.nTraces %#ok<FXUP>
	current = plot(Data.discludedTraces(i,:));
	set(findobj(current,'Type','line'), 'Color', H.traceColor(i,:))
end
hold off

% so there are no problems with subfunctions accessing and updating Data structure
guidata(H.singElecGui, Data)

%% Callbacks for gui

	function moveTracesDummy(hObject, eventdata, H) %#ok<INUSL>
		Data = guidata(H.singElecGui);
		Data = moveTracesManualClustering(hObject, H, Data);
        updateManualClusterAxes(H, Data)
		guidata(H.singElecGui, Data)
    end

    function undoPreviousMove(hObject, eventdata) %#ok<INUSD>
        Data = guidata(H.singElecGui);
        Data.includedIndeces  = Data.prevIncludedIndeces;
        Data.discludedIndeces = Data.prevDiscludedIndeces;
        Data.includedTraces = Data.traces(Data.includedIndeces,:);
        Data.discludedTraces = Data.traces(Data.discludedIndeces,:);
        updateManualClusterAxes(H, Data)
        guidata(H.singElecGui, Data)
	end

	function acceptClusters(hObject, eventdata) %#ok<INUSD>
		Data = guidata(H.singElecGui);
		Data.clusterBin = zeros(1,Data.nTraces);
		Data.clusterBin(Data.includedIndeces) = 1;
		guidata(H.singElecGui, Data)
		close(H.singElecGui)
	end 

uiwait(H.singElecGui)

%sets return value when gui figure is closed by function acceptClusters
clusterBin = Data.clusterBin;

end

