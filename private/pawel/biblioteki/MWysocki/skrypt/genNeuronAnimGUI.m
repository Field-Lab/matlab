function varargout = genNeuronAnimGUI(varargin)
% GENNEURONANIMGUI MATLAB code for genNeuronAnimGUI.fig
%      GENNEURONANIMGUI, by itself, creates a new GENNEURONANIMGUI or raises the existing
%      singleton*.
%
%      H = GENNEURONANIMGUI returns the handle to a new GENNEURONANIMGUI or the handle to
%      the existing singleton*.
%
%      GENNEURONANIMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENNEURONANIMGUI.M with the given input arguments.
%
%      GENNEURONANIMGUI('Property','Value',...) creates a new GENNEURONANIMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before genNeuronAnimGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to genNeuronAnimGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to genNeuronAnimGUI (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help genNeuronAnimGUI

% Last Modified by GUIDE v2.5 01-Dec-2013 02:59:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @genNeuronAnimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @genNeuronAnimGUI_OutputFcn, ...
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


% --- Executes just before genNeuronAnimGUI is made visible.
function genNeuronAnimGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to genNeuronAnimGUI (see VARARGIN)

% Choose default command line output for genNeuronAnimGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes genNeuronAnimGUI wait for user response (see UIRESUME)
% uiwait(handles.window);


% --- Outputs from this function are returned to the command line.
function varargout = genNeuronAnimGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function paramsFile_Callback(hObject, eventdata, handles)
% hObject    handle to paramsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paramsFile as text
%        str2double(get(hObject,'String')) returns contents of paramsFile as a double
    str = get(hObject,'String');
    if ~exist(str,'file')
        set(hObject,'String','<file does not exist>');
        set(hObject,'UserData',[]);
    else
        set(hObject,'UserData',str);
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8 && exist('java.util.HashMap','class') == 8
            file = edu.ucsc.neurobiology.vision.io.ParametersFile(str);
            map = file.getClassIDs();
            tab = unique(cell(map.values().toArray()));
            set(handles.class,'String',tab);
            set(handles.class2,'String',tab);
        end
    end

% --- Executes during object creation, after setting all properties.
function paramsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global param_path;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if ischar(param_path)
        set(hObject,'String',param_path);
        set(hObject,'UserData',param_path);
    else
        set(hObject,'String','');
        set(hObject,'UserData',[]);
    end



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [FileName,PathName] = uigetfile('*.params','Select the Vision parameters file');
    if FileName
        str = [PathName FileName];
        set(handles.paramsFile,'String',str);
        set(handles.paramsFile,'UserData',str);
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8 && exist('java.util.HashMap','class') == 8
            file = edu.ucsc.neurobiology.vision.io.ParametersFile(str);
            map = file.getClassIDs();
            tab = unique(cell(map.values().toArray()));
            set(handles.class,'String',tab);
            set(handles.class2,'String',tab);
        end
    end
    
       


function neuronsFile_Callback(hObject, eventdata, handles)
% hObject    handle to neuronsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuronsFile as text
%        str2double(get(hObject,'String')) returns contents of neuronsFile as a double
    str = get(hObject,'String');
    if ~exist(str,'file')
        set(hObject,'String','<file does not exist>');
        set(hObject,'UserData',[]);
    else
        set(hObject,'UserData',str);
    end

% --- Executes during object creation, after setting all properties.
function neuronsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuronsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global neuron_path;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if ischar(neuron_path)
        set(hObject,'String',neuron_path);
        set(hObject,'UserData',neuron_path);
    else
        set(hObject,'String','');
        set(hObject,'UserData',[]);
    end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [FileName,PathName] = uigetfile('*.neurons','Select the Vision neurons file');
    if FileName
        str = [PathName FileName];
        set(handles.neuronsFile,'String',str);
        set(handles.neuronsFile,'UserData',str);
    end


function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double
    val = str2double(get(hObject,'String'));
    if isnan(val)
        set(hObject,'String','<enter numeric value>');
        set(hObject,'UserData',-1);
    else
        set(hObject,'UserData',val);
    end

% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global time_offset;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(time_offset)
        set(hObject,'String',num2str(time_offset));
        set(hObject,'UserData',time_offset);
    else
        set(hObject,'String','');
        set(hObject,'UserData',-1);
    end



function resolution_Callback(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resolution as text
%        str2double(get(hObject,'String')) returns contents of resolution as a double
    val = str2double(get(hObject,'String'));
    if isnan(val)
        set(hObject,'String','<enter numeric value>');
        set(hObject,'UserData',-1);
    else
        set(hObject,'UserData',val);
    end

% --- Executes during object creation, after setting all properties.
function resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global time_res;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(time_res)
        set(hObject,'String',num2str(time_res));
        set(hObject,'UserData',time_res);
    else
        set(hObject,'String','2.5');
        set(hObject,'UserData',2.5);
    end



function duration_Callback(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of duration as text
%        str2double(get(hObject,'String')) returns contents of duration as a double
    val = str2double(get(hObject,'String'));
    if isnan(val)
        set(hObject,'String','<enter numeric value>');
        set(hObject,'UserData',-1);
    else
        set(hObject,'UserData',val);
    end

% --- Executes during object creation, after setting all properties.
function duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global time_dur;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(time_dur)
        set(hObject,'String',num2str(time_dur));
        set(hObject,'UserData',time_dur);
    else
        set(hObject,'String','1000');
        set(hObject,'UserData',1000);
    end



function class_Callback(hObject, eventdata, handles)
% hObject    handle to class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns class contents as cell array
%        contents{get(hObject,'Value')} returns selected item from class
    contents = cellstr(get(hObject,'String'));
    nclass = contents{get(hObject,'Value')};
    set(hObject,'UserData',nclass);
            
    %pfile = get(handles.paramsFile,'UserData');
    %nclass = get(hObject,'String');
    %if ~isempty(pfile) && ~isempty(get(handles.visionFile,'UserData')) && exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
    %    p = edu.ucsc.neurobiology.vision.io.ParametersFile(pfile);
    %    if isempty(p.getNeuronsInClass(nclass))
    %        set(hObject,'String','<class is not valid>');
    %        set(hObject,'UserData',[]);
    %    else
    %        set(hObject,'UserData',nclass);
    %    end
    %else
    %    set(hObject,'String','<class is unknown>');
    %    set(hObject,'UserData',[]);
    %end
        

% --- Executes during object creation, after setting all properties.
function class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global param_path neuron_classID1;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if ischar(param_path) && exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8 && exist('java.util.HashMap','class') == 8
            file = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);
            map = file.getClassIDs();
            tab = unique(cell(map.values().toArray()));
            set(hObject,'String',tab);
            if ischar(neuron_classID1)
                contents = cellstr(get(hObject,'String'));
                [~,index] = ismember(neuron_classID1,contents);
                set(hObject,'Value',index);
            else
                set(hObject,'Value',1);
            end
    else
        set(hObject,'String','<no classes>');
        set(hObject,'UserData',[]);
    end
    


% --- Executes on button press in double_mode.
function double_mode_Callback(hObject, eventdata, handles)
% hObject    handle to double_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of double_mode
    value = get(hObject,'Value');
    if value == 1
        set(handles.class2,'Enable','on');
        set(hObject,'UserData',1);
    else
        set(handles.class2,'Enable','off');
        set(hObject,'UserData',0);
    end


function class2_Callback(hObject, eventdata, handles)
% hObject    handle to class2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns class2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from class2
    contents = cellstr(get(hObject,'String'));
    nclass = contents{get(hObject,'Value')};
    set(hObject,'UserData',nclass);
    
    %pfile = get(handles.paramsFile,'UserData');
    %nclass = get(hObject,'String');
    %if ~isempty(pfile) && ~isempty(get(handles.visionFile,'UserData')) && exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
    %    p = edu.ucsc.neurobiology.vision.io.ParametersFile(pfile);
    %    if isempty(p.getNeuronsInClass(nclass))
    %        set(hObject,'String','<class is not valid>');
    %        set(hObject,'UserData',[]);
    %    else
    %        set(hObject,'UserData',nclass);
    %    end
    %else
    %    set(hObject,'String','<class is unknown>');
    %    set(hObject,'UserData',[]);
    %end

% --- Executes during object creation, after setting all properties.
function class2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to class2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER. 
    global param_path neuron_classID2 double_mode;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if ischar(param_path) && exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8 && exist('java.util.HashMap','class') == 8
            file = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);
            map = file.getClassIDs();
            tab = unique(cell(map.values().toArray()));
            set(hObject,'String',tab);
            if ischar(neuron_classID2)
                contents = cellstr(get(hObject,'String'));
                [~,index] = ismember(neuron_classID2,contents);
                set(hObject,'Value',index);
            else
                set(hObject,'Value',1);
            end
    else
        set(hObject,'String','<no classes>');
        set(hObject,'UserData',[]);
    end
    if isscalar(double_mode)
        if double_mode == 1
            set(hObject,'Enable','on');
        else
            set(hObject,'Enable','off');
        end
    else
        set(hObject,'Enable','off');
    end



% --- Executes on button press in genNeuronAnimGUI.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to genNeuronAnimGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global generate gif_path gif_width gif_height gif_time_step param_path neuron_path neuron_classID1 neuron_classID2 double_mode time_res time_offset time_dur time_light;
    if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') ~= 8
        w = msgbox('Select Vision library file!','Form is not valid','error');
        return;
    elseif isempty(get(handles.paramsFile,'UserData')) | isempty(get(handles.neuronsFile,'UserData')) | isempty(get(handles.visionFile,'UserData'))
        w = msgbox('Select input file(s)!','Form is not valid','error');
        return;
    elseif get(handles.offset,'UserData') == -1 | get(handles.duration,'UserData') == -1 | get(handles.resolution,'UserData') == -1
        w = msgbox('Specify input time(s)!','Form is not valid','error');
        return;
    elseif isempty(get(handles.class,'UserData')) | (get(handles.double_mode,'UserData') & isempty(get(handles.class2,'UserData')))
        w = msgbox('Specify neurons class(es)!','Form is not valid','error');
        return;
    elseif get(handles.delay,'UserData') == -1 | get(handles.light,'UserData') == -1
        w = msgbox('Specify output time(s)!','Form is not valid','error');
        return;
    else 
        gif_path = get(handles.outputFile,'UserData');
        size = get(handles.gif_size,'UserData');
        gif_width = size(1);
        gif_height = size(2);
        gif_time_step = get(handles.delay,'UserData');
        param_path = get(handles.paramsFile,'UserData');
        neuron_path = get(handles.neuronsFile,'UserData');
        neuron_classID1 = get(handles.class,'UserData');
        neuron_classID2 = get(handles.class2,'UserData');
        double_mode = get(handles.double_mode,'UserData');
        time_res = get(handles.resolution,'UserData'); % [ms] (real time between frames)
        time_offset = get(handles.offset,'UserData'); % [ms]
        time_dur = get(handles.duration,'UserData'); % [ms]
        time_light = get(handles.light,'UserData'); %number of frames when neuron light on
        generate = 1;
        close(handles.window);
    end
    
    
function outputFile_Callback(hObject, eventdata, handles)
% hObject    handle to outputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputFile as text
%        str2double(get(hObject,'String')) returns contents of outputFile as a double
    [pathstr,name,ext] = fileparts(get(hObject,'String'));
    if exist(pathstr,'dir') == 7
        if exist(fullfile(pathstr,name,ext),'dir') == 7
            set(hObject,'String','<selected path is directory>');
            set(hObject,'UserData',[]);
        else
            set(hObject,'UserData',fullfile(pathstr,name));
        end
    else
        set(hObject,'String','<path is not valid>');
        set(hObject,'UserData',[]);
    end


% --- Executes during object creation, after setting all properties.
function outputFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global gif_path;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if ischar(gif_path)
        set(hObject,'String',[gif_path '.gif']);
        set(hObject,'UserData',gif_path);
    else
        set(hObject,'String','');
        set(hObject,'UserData',[]);
    end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [FileName,PathName] = uiputfile('*.gif','Select file for output animation GIF');
    if FileName
        set(handles.outputFile,'String',[PathName FileName]);
        [pathstr,name,ext] = fileparts([PathName FileName]);
        set(handles.outputFile,'UserData',fullfile(pathstr,name));
    end
    

% --- Executes on selection change in gif_size.
function gif_size_Callback(hObject, eventdata, handles)
% hObject    handle to gif_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gif_size contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gif_size
    switch get(hObject,'Value')
        case 1
            set(hObject,'UserData',[1920 1080]);
        case 2
            set(hObject,'UserData',[1280 720]);
        case 3
            set(hObject,'UserData',[1024 576]);
        case 4
            set(hObject,'UserData',[960 540]);
        case 5
            set(hObject,'UserData',[854 480]);
        case 6
            set(hObject,'UserData',[640 360]);
    end


% --- Executes during object creation, after setting all properties.
function gif_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gif_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global gif_width;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(gif_width)
        switch gif_width
            case 1920
                set(hObject,'UserData',[1920 1080]);
                set(hObject,'Value',1);
            case 1280
                set(hObject,'UserData',[1280 720]);
                set(hObject,'Value',2);
            case 1024
                set(hObject,'UserData',[1024 576]);
                set(hObject,'Value',3);
            case 960
                set(hObject,'UserData',[960 540]);
                set(hObject,'Value',4);
            case 854
                set(hObject,'UserData',[854 480]);
                set(hObject,'Value',5);
            case 640
                set(hObject,'UserData',[640 360]);
                set(hObject,'Value',6);
        end
    else
        set(hObject,'UserData',[960 540]);
        set(hObject,'Value',4);
    end



function delay_Callback(hObject, eventdata, handles)
% hObject    handle to delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delay as text
%        str2double(get(hObject,'String')) returns contents of delay as a double
    val = str2double(get(hObject,'String'));
    if isnan(val)
        set(hObject,'String','<enter numeric value>');
        set(hObject,'UserData',-1);
    else
        if val < 0.06
            val = 0.06;
            set(hObject,'String','0.06');
        end
        set(hObject,'UserData',val);
    end

% --- Executes during object creation, after setting all properties.
function delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global gif_time_step;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(gif_time_step)
        set(hObject,'String',num2str(gif_time_step));
        set(hObject,'UserData',gif_time_step);
    else
        set(hObject,'String','0.06');
        set(hObject,'UserData',0.06);
    end


function light_Callback(hObject, eventdata, handles)
% hObject    handle to light (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of light as text
%        str2double(get(hObject,'String')) returns contents of light as a double
    val = str2num(get(hObject,'String'));
    if isempty(val)
        set(hObject,'String','<enter numeric value>');
        set(hObject,'UserData',-1);
    else
        val(val < 1) = 1;
        rval = round(val(1));
        set(hObject,'String',num2str(rval));
        set(hObject,'UserData',rval);
    end

% --- Executes during object creation, after setting all properties.
function light_CreateFcn(hObject, eventdata, handles)
% hObject    handle to light (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    global time_light;
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if isscalar(time_light)
        set(hObject,'String',num2str(time_light));
        set(hObject,'UserData',time_light);
    else
        set(hObject,'String','2');
        set(hObject,'UserData',2);
    end



function visionFile_Callback(hObject, eventdata, handles)
% hObject    handle to visionFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visionFile as text
%        str2double(get(hObject,'String')) returns contents of visionFile as a double
    str = get(hObject,'String');
    if ~exist(str,'file')
        set(hObject,'String','<file does not exist>');
        set(hObject,'UserData',[]);
    else
        set(hObject,'UserData',str);
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') ~= 8
            javaaddpath(str);
        end
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
            set(hObject,'Enable','off');
            set(handles.loaded,'Value',1);
            set(handles.loaded,'Enable','off');
            set(handles.visionButton,'Enable','off');
            
            pfile = get(handles.paramsFile,'UserData');
            if ~isempty(pfile) && exist('java.util.HashMap','class') == 8
                file = edu.ucsc.neurobiology.vision.io.ParametersFile(pfile);
                map = file.getClassIDs();
                tab = unique(cell(map.values().toArray()));
                set(handles.class,'String',tab);
                set(handles.class2,'String',tab);
            end
        end
    end

% --- Executes during object creation, after setting all properties.
function visionFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visionFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
        set(hObject,'Enable','off');
        set(hObject,'UserData','loaded');
    else
        set(hObject,'Enable','on');
        set(hObject,'UserData',[]);
    end


% --- Executes on button press in visionButton.
function visionButton_Callback(hObject, eventdata, handles)
% hObject    handle to visionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [FileName,PathName] = uigetfile('*.jar','Select the Vision jar library');
    if FileName
        str = [PathName FileName];
        set(handles.visionFile,'String',str);
        set(handles.visionFile,'UserData',str);
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') ~= 8
            javaaddpath(str);
        end
        if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
            set(handles.visionFile,'Enable','off');
            set(handles.loaded,'Value',1);
            set(handles.loaded,'Enable','off');
            set(hObject,'Enable','off');
            
            pfile = get(handles.paramsFile,'UserData');
            if ~isempty(pfile) && exist('java.util.HashMap','class') == 8
                file = edu.ucsc.neurobiology.vision.io.ParametersFile(pfile);
                map = file.getClassIDs();
                tab = unique(cell(map.values().toArray()));
                set(handles.class,'String',tab);
                set(handles.class2,'String',tab);
            end
        end
    end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in loaded.
function loaded_Callback(hObject, eventdata, handles)
% hObject    handle to loaded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loaded
    if get(hObject,'Value')
        set(handles.visionFile,'Enable','off');
        set(handles.visionFile,'UserData','loaded');
        set(handles.visionButton,'Enable','off');
    else
        set(handles.visionFile,'Enable','on');
        set(handles.visionButton,'Enable','on');
    end


% --- Executes during object creation, after setting all properties.
function loaded_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loaded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
        set(hObject,'Value',1);
        set(hObject,'Enable','off');
    else
        set(hObject,'Value',0);
        set(hObject,'Enable','on');
    end


% --- Executes during object creation, after setting all properties.
function visionButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    if exist('edu.ucsc.neurobiology.vision.io.ParametersFile','class') == 8
        set(hObject,'Enable','off');
    else
        set(hObject,'Enable','on');
    end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    m = msgbox('This program was written by Mateusz Wysocki, the AGH UST student, as part of his Engineering Thesis in 2013.','About');


% --- Executes during object creation, after setting all properties.
function double_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to double_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    global double_mode;
    if isscalar(double_mode)
        if double_mode == 1
            set(hObject,'Value',1);
            set(hObject,'UserData',1);
        else
            set(hObject,'Value',0);
            set(hObject,'UserData',0);
        end
    else
        set(hObject,'Value',0);
        set(hObject,'UserData',0);
    end


% --- Executes during object creation, after setting all properties.
function window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    global generate;
    generate = 0;


% --- Executes when user attempts to close window.
function window_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
