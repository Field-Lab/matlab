function varargout = WaveformClassifier1(varargin)
% WAVEFORMCLASSIFIER1 M-file for WaveformClassifier1.fig
%      WAVEFORMCLASSIFIER1, by itself, creates a new WAVEFORMCLASSIFIER1 or raises the existing
%      singleton*.
%
%      H = WAVEFORMCLASSIFIER1 returns the handle to a new WAVEFORMCLASSIFIER1 or the handle to
%      the existing singleton*.
%
%      WAVEFORMCLASSIFIER1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVEFORMCLASSIFIER1.M with the given input arguments.
%
%      WAVEFORMCLASSIFIER1('Property','Value',...) creates a new WAVEFORMCLASSIFIER1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WaveformClassifier1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WaveformClassifier1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WaveformClassifier1

% Last Modified by GUIDE v2.5 10-Mar-2014 15:18:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WaveformClassifier1_OpeningFcn, ...
                   'gui_OutputFcn',  @WaveformClassifier1_OutputFcn, ...
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


% --- Executes just before WaveformClassifier1 is made visible.
function WaveformClassifier1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WaveformClassifier1 (see VARARGIN)

% Choose default command line output for WaveformClassifier1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WaveformClassifier1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WaveformClassifier1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [NewEventID,NewNeuronID,NewElectrodeID]=ModifyFile(EventID,Value,handles)%(eventdata,EventID,Value,hObject,handles);
Data=ReadData(handles);
Data(3,EventID)=Value
FolderPath=get(handles.edit1,'String');
fid = fopen([FolderPath '\NeuronsPatterns2.bin'],'wb');
fwrite(fid,Data,'int32');
fclose(fid);
NewEventID=EventID+1;
NewNeuronID=Data(1,NewEventID);
NewElectrodeID=Data(2,NewEventID);

function r=ShowFigure(NeuronID,ElectrodeID,handles);
FolderPath=get(handles.edit1,'String');
FigureName=[FolderPath '\Neuron' num2str(NeuronID) '_Stim' num2str(ElectrodeID) 'add.tif']
RGB = imread(FigureName);
[m,n]=size(RGB);
figure(1)
set(gcf,'Units','pixels','Position',[0 0 n/3*1.111 m*1.111])
image(RGB); %colormap(map)
set(gca,'Position',[0.05 0.05 0.9 0.9])
set(gcf,'NextPlot','replace');
r=2;

function EventID=NextEvent(Value,handles);
EventID=str2num(get(handles.edit2,'String'))
[NewEventID,NewNeuronID,NewElectrodeID]=ModifyFile(EventID,Value,handles)
r=ShowFigure(NewNeuronID,NewElectrodeID,handles);
set(handles.edit2,'String',num2str(NewEventID));
set(handles.edit3,'String',num2str(NewNeuronID));
set(handles.edit4,'String',num2str(NewElectrodeID));

function Data=ReadData(handles);
FolderPath=get(handles.edit1,'String');
fid = fopen([FolderPath '\NeuronsPatterns2.bin'],'r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
EventID=NextEvent(1,handles);
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
EventID=NextEvent(2,handles);

% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
EventID=NextEvent(3,handles);
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
EventID=NextEvent(4,handles);
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
Data=ReadData(handles);
set(handles.edit2,'String',num2str(1));
set(handles.edit3,'String',Data(1,1));
set(handles.edit4,'String',Data(2,1));
r=ShowFigure(Data(1,1),Data(2,1),handles);

% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
EventID=str2num(get(handles.edit2,'String'))
Data=ReadData(handles);
NeuronID=Data(1,EventID);
ElectrodeID=Data(2,EventID);
set(handles.edit3,'String',NeuronID);
set(handles.edit4,'String',ElectrodeID);
r=ShowFigure(NeuronID,ElectrodeID,handles);
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
