function varargout = StimScanWatch5b(varargin)
% STIMSCANWATCH5B M-file for StimScanWatch5b.fig
%      STIMSCANWATCH5B, by itself, creates a new STIMSCANWATCH5B or raises
%      the existing
%      singleton*.
%
%      H = STIMSCANWATCH5B returns the handle to a new STIMSCANWATCH5B or the
%      handle to
%      the existing singleton*.
%
%      STIMSCANWATCH5B('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in STIMSCANWATCH5B.M with the given input
%      arguments.
%
%      STIMSCANWATCH5B('Property','Value',...) creates a new STIMSCANWATCH5B
%      or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before StimScanWatch5b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to StimScanWatch5b_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StimScanWatch5b

% Last Modified by GUIDE v2.5 26-Sep-2012 17:03:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StimScanWatch5b_OpeningFcn, ...
                   'gui_OutputFcn',  @StimScanWatch5b_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before StimScanWatch5b is made visible.
function StimScanWatch5b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StimScanWatch5b (see VARARGIN)

% Choose default command line output for StimScanWatch5b
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StimScanWatch5b wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = StimScanWatch5b_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% --- Executes on selection change in popupmenu1.

function edit1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

function edit1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton4_Callback(hObject, eventdata, handles)
%"Next movie" button
MovieNumber=str2num(get(handles.edit1,'String'));
step=str2num(get(handles.edit7,'String'));
set(handles.edit1,'String',num2str(MovieNumber+step));
[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);
ShowArtifacts=get(handles.radiobutton2,'Value');
UpdateTraces=get(handles.radiobutton10,'Value');
UpdateEI=get(handles.radiobutton11,'Value');
if ShowArtifacts==1
    Traces=cat(1,DataTraces,ArtifactDataTraces);
    SADT=size(ArtifactDataTraces);
    WaveformTypesFinal=cat(1,WaveformTypes,zeros(SADT(1),1));
else
    Traces=DataTraces;
    WaveformTypesFinal=WaveformTypes;
end
if UpdateTraces==1
       y=PlotData(Traces,Channels,WaveformTypesFinal,hObject, eventdata, handles);
end
if UpdateEI==1
        y=PlotEI(Traces,Channels,WaveformTypes,hObject, eventdata, handles);
end
%y=PlotData(Traces,Channels,WaveformTypesFinal,hObject, eventdata, handles);

function edit2_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a
%        double


function edit2_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1

function edit6_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton5_Callback(hObject, eventdata, handles)
%Previous Movie button
MovieNumber=str2num(get(handles.edit1,'String'));
step=str2num(get(handles.edit7,'String'));
NextMovieNumber=MovieNumber-step;
UpdateTraces=get(handles.radiobutton10,'Value');
UpdateEI=get(handles.radiobutton11,'Value');
if NextMovieNumber>0
    set(handles.edit1,'String',num2str(MovieNumber-step));
    MovieNumber=str2num(get(handles.edit1,'String'));
    [DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
    WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);
    
    ShowArtifacts=get(handles.radiobutton2,'Value');
    if ShowArtifacts==1
        Traces=cat(1,DataTraces,ArtifactDataTraces);
        SADT=size(ArtifactDataTraces);
        WaveformTypesFinal=cat(1,WaveformTypes,zeros(SADT(1),1));
    else
        Traces=DataTraces;
        WaveformTypesFinal=WaveformTypes;
    end
    if UpdateTraces==1
        y=PlotData(Traces,Channels,WaveformTypesFinal,hObject, eventdata, handles);
    end
    if UpdateEI==1
        y=PlotEI(Traces,Channels,WaveformTypes,hObject, eventdata, handles);
    end
    %y=PlotData(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles);
end

function pushbutton6_Callback(hObject, eventdata, handles)
%Refresh button

[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);

ShowArtifacts=get(handles.radiobutton2,'Value');
if ShowArtifacts==1
    Traces=cat(1,DataTraces,ArtifactDataTraces);
    SADT=size(ArtifactDataTraces);
    WaveformTypesFinal=cat(1,WaveformTypes,zeros(SADT(1),1));
else
    Traces=DataTraces;
    WaveformTypesFinal=WaveformTypes;
end
y=PlotData(Traces,Channels,WaveformTypesFinal,hObject, eventdata, handles);

function Channels=GUI_Channels(handles)
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
IfPattern=get(handles.radiobutton6,'Value');
if (IfPattern==1)
    PatternNumber=str2num(get(handles.edit2,'String'));
    CenterChannel=PatternNumber;
else
    CenterChannel=str2num(get(handles.edit14,'String'));
end
Radius=str2num(get(handles.edit15,'String'))
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
b=find(Channels~=9 & Channels~=25 & Channels~=57 & Channels~=4);
Channels=Channels(b);

function pushbutton7_Callback(hObject, eventdata, handles)
%"Cluster traces" button
[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
MovieNumber=str2num(get(handles.edit1,'String'));
PatternNumber=str2num(get(handles.edit2,'String'));
SPfilename=get(handles.edit13,'String');
Pattern=ReadPattern(hObject, eventdata, handles);

[ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
IfPattern=get(handles.radiobutton6,'Value');
if (IfPattern==1)
    CenterChannel=PatternNumber;
else
    CenterChannel=str2num(get(handles.edit14,'String'));
end
Radius=str2num(get(handles.edit15,'String'));
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
OutputTraces=NS_CreateVectorForPCA2(DataTraces,Channels,Pattern,NS_GlobalConstants);
figure(2)
plot(OutputTraces');
% remove all the samples when given channel was disconnected!!
WaveformTypes = manualPCACluster(OutputTraces);
size(WaveformTypes)

ShowArtifacts=get(handles.radiobutton2,'Value');
if ShowArtifacts==1
    Traces=cat(1,DataTraces,ArtifactDataTraces);
    SADT=size(ArtifactDataTraces)
    WaveformTypesFinal=cat(1,WaveformTypes,zeros(SADT(1),1));
else
    Traces=DataTraces;
    WaveformTypesFinal=WaveformTypes;
end
%Channels=[16:19];
y=PlotData(Traces,Channels,WaveformTypesFinal,hObject, eventdata, handles);
%y=PlotData(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles);
ClusterFilename=get(handles.edit10,'String');

SaveWaveformTypes=get(handles.radiobutton4,'Value');
if SaveWaveformTypes==1
    PatternNumber=str2num(get(handles.edit2,'String'));
    MovieNumber=str2num(get(handles.edit1,'String'));
    FilePath=get(handles.edit10,'String');
    header=NS_SaveClusterFile(FilePath,PatternNumber,MovieNumber,WaveformTypes);
    %WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber);
end

function pushbutton8_Callback(hObject, eventdata, handles);
%Show stimuli
SPfilename=get(handles.edit13,'String');
MovieNumber=str2num(get(handles.edit1,'String'))
PatternNumber=str2num(get(handles.edit2,'String'))

[ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles);
NS_GlobalConstants=NS_GenerateGlobalConstants(ArrayCount);
%ArrayID=1;

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
IfPattern=get(handles.radiobutton6,'Value');
if (IfPattern==1)
    CenterChannel=PatternNumber;
else
    CenterChannel=str2num(get(handles.edit14,'String'));
end
Radius=str2num(get(handles.edit15,'String'));
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
%channels=electrodeMap.getAdjacentsTo(patternNumberToPlot,1)';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-0.6 0.6],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','current [\muA]');

Status=ReadStatus(hObject, eventdata, handles);
Pattern=ReadPattern(hObject, eventdata, handles);
Channels
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
figHandle = NS_PlotClustersOfSignaturesOnArrayLayout(traces,Channels,0,ArrayID,FigureProperties,NS_GlobalConstants,patternTimes);
%figHandle = plotStimulusTraces(patternNumberToPlot, channels, ArrayID, FigureProperties, NS_GlobalConstants);

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
% str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Status=ReadStatus(hObject, eventdata, handles)
DataPath=get(handles.edit3,'String');
%PatternNumber=str2num(get(handles.edit2,'String'));
MovieNumber=str2num(get(handles.edit1,'String'));
FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
%clear Pattern;
load(FullName_status);

function Pattern=ReadPattern(hObject, eventdata, handles);
DataPath=get(handles.edit3,'String');
PatternNumber=str2num(get(handles.edit2,'String'));
MovieNumber=str2num(get(handles.edit1,'String'));
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
clear Pattern;
load(FullName_pattern);
size(Pattern)

function [ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles);
a=get(handles.popupmenu2,'Value');

switch a
    case 1
        ArrayID=1;
        ArrayCount=61;
    case 2
        ArrayID=500;
        ArrayCount=512;
    otherwise
        ArrayID=1502;
        ArrayCount=519;
end

function [DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
%Reads data from small files
[ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles)

ArtifactSubtraction=get(handles.radiobutton1,'Value');
NS_GlobalConstants=NS_GenerateGlobalConstants(ArrayCount);

DataPath=get(handles.edit3,'String');
ArtifactDataPath=get(handles.edit6,'String');
PatternNumber=str2num(get(handles.edit2,'String'));
MovieNumber=str2num(get(handles.edit1,'String'));

LowLimit=str2num(get(handles.edit8,'String'));
HighLimit=str2num(get(handles.edit9,'String'));

IfAllEvents=get(handles.radiobutton12,'Value');
if IfAllEvents==1
    EventNumber=0;
else
    EventNumber0=str2num(get(handles.edit23,'String'));
    EventNumber1=str2num(get(handles.edit31,'String'));
    EventNumber=[EventNumber0:EventNumber1]
end
DataPath
PatternNumber
MovieNumber
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,1500,EventNumber);
rrrr=size(DataTraces)

IfPrevious=get(handles.radiobutton13,'Value');
if (IfPrevious==1)
    [DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber-1,100,EventNumber);
end

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
IfPattern=get(handles.radiobutton6,'Value');
if (IfPattern==1)
    CenterChannel=PatternNumber;
else
    CenterChannel=str2num(get(handles.edit14,'String'));
end
Radius=str2num(get(handles.edit15,'String'));
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)'
%Channels=[16:19];
if ArrayCount==61
    b=find(Channels~=9 & Channels~=25 & Channels~=57 & Channels~=4);
    Channels=Channels(b);
end
TraceLength=str2num(get(handles.edit30,'String'));
if (IfPrevious==1)
    DataTracesFin=zeros(80,length(Channels),TraceLength);
    DataTracesFin(1:40,:,:)=DataTraces(1:40,Channels,1:TraceLength);
    DataTracesFin(41:80,:,:)=DataTraces2(1:40,Channels,1:TraceLength);
    DataTraces=DataTracesFin;    
else
    size(DataTraces)
    Channels
    TraceLength
    DataTraces=DataTraces(:,Channels,1:TraceLength);
end
ArtifactDataTraces=ArtifactDataTraces(:,Channels,:);

function WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);
ReadClusterIndex=get(handles.radiobutton3,'Value');
SD=size(DataTraces);
if ReadClusterIndex==1
    PatternNumber=str2num(get(handles.edit2,'String'));
    MovieNumber=str2num(get(handles.edit1,'String'));
    FileName=get(handles.edit10,'String');
    WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber,SD(1));
else
    %SD=size(DataTraces);
    WaveformTypes=ones(SD(1),1);
end

function y=PlotData(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles)
size(DataTraces)
[ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles);
LowLimit=str2num(get(handles.edit8,'String'));
HighLimit=str2num(get(handles.edit9,'String'));
NumberOfSamples=str2num(get(handles.edit30,'String'));
FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 NumberOfSamples],'AmplitudeRange',[LowLimit HighLimit],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
figure(FigureProperties.FigureNumber);
clf;
NS_GlobalConstants=NS_GenerateGlobalConstants(ArrayCount);
SubtractCluster=get(handles.radiobutton7,'Value');
if SubtractCluster==1
    ClusterToSubtract=str2num(get(handles.edit16,'String'));
    a=find(WaveformTypes==ClusterToSubtract);
    TraceToSubtract=mean(DataTraces(a,:,:));
    ST=size(DataTraces);
    for i=1:ST(1)
        DataTraces(i,:,:)=DataTraces(i,:,:)-TraceToSubtract;
    end
end

IfPrevious=get(handles.radiobutton13,'Value');
if IfPrevious==1
    WaveformTypes(1:40)=1;
    WaveformTypes(41:80)=2;
end
ClusterVsEventPlot=get(handles.radiobutton14,'Value');
if ClusterVsEventPlot==1
    figure(111);
    h=plot(WaveformTypes,'bd-');
    set(h,'MarkerSize',10);
    h=gca;
    set(h,'FontSize',24);
    set(h,'YLim',[0.8 2.2])
    set(h,'YTick',[1 2])
    set(h,'YTickLabel',{'false' 'success'});
    set(h,'XTick',[0:10:100]);
    hxl=xlabel('trial index');
    set(hxl,'FontSize',24);
    grid on;
end
khjkhjk=size(DataTraces)
y=NS_PlotClustersOfSignaturesOnArrayLayout(DataTraces,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);

function y=PlotEI(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles);
a=get(handles.popupmenu2,'Value');

switch a
    case 1
        NS_GlobalConstants=NS_GenerateGlobalConstants(61);
    case 2
        NS_GlobalConstants=NS_GenerateGlobalConstants(512);
    otherwise
        NS_GlobalConstants=NS_GenerateGlobalConstants(512);
end

LowLimit=str2num(get(handles.edit24,'String'));
HighLimit=str2num(get(handles.edit25,'String'));
FigureProperties=struct('FigureNumber',4,'Subplot',[2 3 3],'TimeRange',[0 100],'AmplitudeRange',[LowLimit HighLimit],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','signal [mV]');
%NS_GlobalConstants=NS_GenerateGlobalConstants(61);

ClusterToSubtract=str2num(get(handles.edit16,'String'));
a=find(WaveformTypes==ClusterToSubtract);
TraceToSubtract=mean(DataTraces(a,:,:));
    
ClusterToShow=str2num(get(handles.edit17,'String'));
a=find(WaveformTypes==ClusterToShow);
TracesToShow=mean(DataTraces(a,:,:))-TraceToSubtract;
%TracesToShow=mean(DataTraces(:,:,:)) 

WaveformTypes=[1];
ShowNeuronEI=get(handles.radiobutton8,'Value');
if ShowNeuronEI==1    
    ShowNeuronID=str2num(get(handles.edit18,'String'));
    ST=size(TracesToShow);
    NeuronID=str2num(get(handles.edit18,'String'));
    EI=NeuronSignature(NeuronID,ST(3),handles);
    [EIshift,MaxCorr]=NS_FindShiftBetweenEIs(reshape(TracesToShow(1,:,:),ST(2),ST(3)),EI);
    EI2=EI(:,EIshift+1:EIshift+ST(3));        
    TracesToShow(2,1:ST(2),1:ST(3))=EI2;
    WaveformTypes=[1 2];
end
[ArrayID,ArrayCount]=ReadChannelsCount(hObject, eventdata, handles);
y=NS_PlotClustersOfSignaturesOnArrayLayout(TracesToShow,Channels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);

function EI=NeuronEI(NeuronID,TracesLength,handles);
Channels=GUI_Channels(handles);
RAWFileName=get(handles.edit19,'String');
%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(get(handles.edit20,'String'));
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(get(handles.edit21,'String'));
eiFile=edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(get(handles.edit20,'String'));
EI0=eiFile.getImage(NeuronID);
size(EI0)
EI=reshape(EI0(1,Channels+1,1:50),length(Channels),50);


function EI=NeuronSignature(NeuronID,TracesLength,handles);
%The function finds EI and subtracts offest
a=get(handles.popupmenu2,'Value');

switch a
    case 1
        NS_GlobalConstants=NS_GenerateGlobalConstants(61);
    case 2
        NS_GlobalConstants=NS_GenerateGlobalConstants(512);
    otherwise
        NS_GlobalConstants=NS_GenerateGlobalConstants(512);
end

Channels=GUI_Channels(handles);
RAWFileName=get(handles.edit19,'String');
paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(get(handles.edit20,'String'));
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(get(handles.edit21,'String'));
%idList = neuronFile.getIDList();

spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
Timings=spikeTimes(1,1:min(100,length(spikeTimes)))-20;
TracesLength=100;
[RAWtraces,EI]=NS_AverageTraces(RAWFileName,Timings,Channels,[1 TracesLength],NS_GlobalConstants);
SEI=size(EI);
for i=1:SEI(2)
    EI(:,i)=EI(:,i)-(EI(1,i)+EI(SEI(1),i))/2;
end
EI=EI';

function NeuronInfo=ReadNeuronInfoFromThresholdFile(TresholdFilePath,NeuronID,handles);
fid1=fopen(TresholdFilePath,'r')
b=fread(fid1,inf,'double')
fclose(fid1)
c11=reshape(b,length(b)/4,4)
NeuronInfoCoordinate=find(c11(:,1)==NeuronID)
NeuronInfo=c11(NeuronInfoCoordinate,:);

function SaveNeuronInfoIntoThresholdFile(TresholdFilePath,NeuronID,PrimaryStimulatingElectrode,ThresholdCurrent);
fid1=fopen(TresholdFilePath,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
c11=reshape(b,length(b)/4,4);
NeuronInfoCoordinate=find(c11(:,1)==NeuronID);
NeuronInfo=c11(NeuronInfoCoordinate,:)
NeuronInfo(1,3)=PrimaryStimulatingElectrode;
NeuronInfo(1,4)=ThresholdCurrent;
c11(NeuronInfoCoordinate,:)=NeuronInfo;
fid1=fopen(TresholdFilePath,'w');
b=fwrite(fid1,c11,'double');
fclose(fid1);

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
%function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Traces=ReadData(hObject, eventdata, handles)
%clusterIndex = manualPCACluster(Traces)
%refresh(hObject, eventdata, handles);



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton10_Callback(hObject, eventdata, handles)
% Show EI
[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);

ShowArtifacts=get(handles.radiobutton2,'Value');
if ShowArtifacts==1
    Traces=cat(1,DataTraces,ArtifactDataTraces);
    SADT=size(ArtifactDataTraces);
    WaveformTypesFinal=cat(1,WaveformTypes,zeros(SADT(1),1));
else
    Traces=DataTraces;
    WaveformTypesFinal=WaveformTypes;
end
save Traces Traces;
y=PlotEI(Traces,Channels,WaveformTypes,hObject, eventdata, handles);

function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a
%        double

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a
%        double

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
%Match Neuron
%1. Find differantial EI of stimulated neuron
[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
WaveformTypes=ImportWaveformTypes(DataTraces,hObject, eventdata, handles);

ClusterToSubtract=str2num(get(handles.edit16,'String'));
a=find(WaveformTypes==ClusterToSubtract);
TraceToSubtract=mean(DataTraces(a,:,:));
    
ClusterToShow=str2num(get(handles.edit17,'String'));
a=find(WaveformTypes==ClusterToShow);
TracesToShow=mean(DataTraces(a,:,:))-TraceToSubtract; %this is the differential EI
ST=size(TracesToShow);
    
Channels=GUI_Channels(handles);
RAWFileName=get(handles.edit19,'String');

eiFile=edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(get(handles.edit20,'String'));
%eiFile=edu.ucsc.neurobiology.vision.io.ParametersFile(get(handles.edit20,'String'));
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(get(handles.edit21,'String'));
idList = neuronFile.getIDList();
MaxCorrs=zeros(1,length(idList));

TracesLength=100;
for i=1:length(idList)
    NeuronID=idList(i)
    %EI=NeuronSignature(NeuronID,TracesLength,handles);
    %eis=size(EI)
    EI=NeuronEI(NeuronID,TracesLength,handles);
    eis=size(EI)
    NormEI=NS_NormalizeEI(EI);   
    
    %EI0=eiFile.getImage(NeuronID);
    %size(EI0)
    %EI=reshape(EI0(1,TracesToShow,:),length(TracesToShow),81);
    
    [EIshift,MaxCorr]=NS_FindShiftBetweenEIs(reshape(TracesToShow(1,:,:),ST(2),ST(3)),NormEI);
    MaxCorrs(1,i)=MaxCorr;
end
figure(122)
plot(idList,MaxCorrs,'bd-');
grid on;
M=max(MaxCorrs)
IDindex=find(MaxCorrs==M)
set(handles.edit18,'String',num2str(idList(IDindex)));


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton11



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton13



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(handles.popupmenu2,'Value');
set(handles.text24,'String',a);


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[DataTraces,ArtifactDataTraces,Channels,NS_GlobalConstants]=ReadData(hObject, eventdata, handles);
SDT=size(DataTraces);
MovieNumber=str2num(get(handles.edit1,'String'));
PatternNumber=str2num(get(handles.edit2,'String'));
WaveformTypes=ones(1,SDT(1))*2;
y=PlotData(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles);
%y=PlotData(DataTraces,Channels,WaveformTypes,hObject, eventdata, handles);
ClusterFilename=get(handles.edit10,'String');

SaveWaveformTypes=get(handles.radiobutton4,'Value');
if SaveWaveformTypes==1
    PatternNumber=str2num(get(handles.edit2,'String'));
    MovieNumber=str2num(get(handles.edit1,'String'));
    FilePath=get(handles.edit10,'String');
    header=NS_SaveClusterFile(FilePath,PatternNumber,MovieNumber,WaveformTypes);
    %WaveformTypes=NS_ReadClusterFile(FileName,MovieNumber,PatternNumber);
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TresholdFilePath=get(handles.edit13,'String');
PrimaryStimulatingElectrode=str2num(get(handles.edit2,'String'));
MovieNumber=str2num(get(handles.edit1,'String'));
DataPath=get(handles.edit3,'String');
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,PrimaryStimulatingElectrode,PrimaryStimulatingElectrode,MovieNumber,NS_GlobalConstants)
NeuronID=str2num(get(handles.edit34,'String'));
SaveNeuronInfoIntoThresholdFile(TresholdFilePath,NeuronID,PrimaryStimulatingElectrode,Amplitude);

function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles);
TresholdFilePath=get(handles.edit13,'String');
NeuronID=str2num(get(handles.edit34,'String'));
NeuronInfo=ReadNeuronInfoFromThresholdFile(TresholdFilePath,NeuronID,handles)
set(handles.edit14,'String',num2str(NeuronInfo(1,2)));
PrimaryStim=NeuronInfo(1,3);
if PrimaryStim~=0
    set(handles.edit2,'String',num2str(NeuronInfo(1,3)));
else
    set(handles.edit2,'String',num2str(NeuronInfo(1,2)));
end
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
