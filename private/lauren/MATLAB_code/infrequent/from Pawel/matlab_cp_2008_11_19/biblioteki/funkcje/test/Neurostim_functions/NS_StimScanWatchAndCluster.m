function h4=NS_StimScanWatchAndCluster;

LeftMarg=50;
TopMarg=50;

XStep=60;
YStep=20;

WindowWidth=600;
WindowHeight=300;
h.gui = figure('Visible', 'on', 'Position',[100,100,WindowWidth,WindowHeight]);
h1=uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'Next Movie', 'Position', [LeftMarg WindowHeight-TopMarg-YStep*9 100 30], 'Callback', @PC1);
%get(h1)
h2=uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'Previous Movie', 'Position', [LeftMarg+120 WindowHeight-TopMarg-YStep*9 100 30], 'Callback', @PC2);

h3=uicontrol(h.gui, 'Style', 'text', 'String', 'Data Path:', 'Position', [LeftMarg WindowHeight-TopMarg 100 14]);
h4=uicontrol(h.gui, 'Style', 'edit', 'Position', [LeftMarg+110 WindowHeight-TopMarg 250 14]);
get(h4)
set(h4,'string','ergerery');
guidata(h.gui,h4);

uicontrol(h.gui, 'Style', 'text', 'String', 'Artifact Data Path:', 'Position', [LeftMarg WindowHeight-TopMarg-YStep  100 14]);
uicontrol(h.gui, 'Style', 'edit', 'Position', [LeftMarg+110 WindowHeight-TopMarg-YStep 250 14]);
uicontrol(h.gui, 'Style', 'radiobutton', 'String', 'Subtract Artifacts', 'Position', [LeftMarg+370 WindowHeight-TopMarg-YStep 120 14]);

uicontrol(h.gui, 'Style', 'text', 'String', 'Movie Number:', 'Position', [LeftMarg WindowHeight-TopMarg-YStep*3 100 14]);
uicontrol(h.gui, 'Style', 'edit', 'Position', [LeftMarg+110 WindowHeight-TopMarg-YStep*3 30 14]);

uicontrol(h.gui, 'Style', 'text', 'String', 'PatternNumber:', 'Position', [LeftMarg WindowHeight-TopMarg-YStep*4 100 14]);
uicontrol(h.gui, 'Style', 'edit', 'Position', [LeftMarg+110 WindowHeight-TopMarg-YStep*4 30 14]);

uicontrol(h.gui, 'Style', 'text', 'String', 'Movie Index Step:', 'Position', [LeftMarg WindowHeight-TopMarg-YStep*6 100 14]);
uicontrol(h.gui, 'Style', 'edit', 'Position', [LeftMarg+110 WindowHeight-TopMarg-YStep*6 30 14]);

data=guihandles(h.gui)
y=1;
%end

function PC1(hObject, eventdata)
    'aaa'
    'bbb'
    
%end

NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArrayID=1;

FigureProperties=struct('FigureNumber',2,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-1200 450],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
Patterns=[1];
PatternsIndexes=[2];

get(h4)
DataPath=get(h4,'String')
ArtifactDataPath=get(handles.edit2,'String')
PatternNumber=str2num(get(handles.edit4,'String'))
MovieNumber=str2num(get(handles.edit3,'String'))

figure(FigureProperties.FigureNumber);
clf;
figure(FigureProperties.FigureNumber+10);
clf;
%cd(DataPath);

[Channels]=NS_ShowPreprocessedData(DataPath,ArtifactDataPath,PatternNumber,MovieNumber,Patterns,PatternsIndexes,FigureProperties,NS_GlobalConstants);
%end

function PC2(hObject, eventdata)
close
%end