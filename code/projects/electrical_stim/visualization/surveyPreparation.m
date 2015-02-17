function varargout = surveyPreparation(varargin)
% SURVEYPREPARATION MATLAB code for surveyPreparation.fig
%      SURVEYPREPARATION, by itself, creates a new SURVEYPREPARATION or raises the existing
%      singleton*.
%
%      H = SURVEYPREPARATION returns the handle to a new SURVEYPREPARATION or the handle to
%      the existing singleton*.
%
%      SURVEYPREPARATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SURVEYPREPARATION.M with the given input arguments.
%
%      SURVEYPREPARATION('Property','Value',...) creates a new SURVEYPREPARATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before surveyPreparation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to surveyPreparation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help surveyPreparation

% Last Modified by GUIDE v2.5 12-Feb-2015 14:53:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @surveyPreparation_OpeningFcn, ...
                   'gui_OutputFcn',  @surveyPreparation_OutputFcn, ...
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


% --- Executes just before surveyPreparation is made visible.
function surveyPreparation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to surveyPreparation (see VARARGIN)

% Choose default command line output for surveyPreparation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes surveyPreparation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = surveyPreparation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable1,'Data',[]);
handles.dataPath = get(handles.datapath,'String'); 
datarun  = load_data(handles.dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, 'all');
handles.datarun = datarun; 

handles.spikeThresh = str2double(get(handles.threshold, 'String'));
[cellIdsToCheck, cellIndices, electrodes] = getLargeAmpSpikes(handles.datarun,handles.spikeThresh); 
% 3. Determine firing rate stability of these cells
N = length(cellIdsToCheck); 
T = datarun.duration; % length of run

binsize     =  5; % seconds
binranges   = binsize:binsize:T; 
firing_rate = zeros(N,1);
firing_rate_vector = zeros(N,length(binranges)); 

for i=1:N
    firing_rate(i) = length(datarun.spikes{cellIndices(i)})/T;
    bincounts = histc(datarun.spikes{cellIndices(i)},binranges); 
    firing_rate_vector(i,:) = bincounts/binsize; 
end
devs                = std(firing_rate_vector,[],2);
handles.binranges   = binranges; 
handles.firingRates = firing_rate_vector; 

goodCol = num2cell(false(N,1));
data = num2cell(cat(2,cellIdsToCheck',electrodes,firing_rate,devs));
% toCat = num2cell(goodCol); 
% data = {(cellIdsToCheck'),(electrodes),(firing_rate),(devs),(goodCol)};
% data = {num2cell(cellIdsToCheck'),num2cell(electrodes),num2cell(firing_rate),num2cell(devs),num2cell(goodCol)};
data = [data goodCol]; 
set(handles.uitable1,'Data',data);
% Update handles structure
guidata(hObject, handles);


function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double
update_Callback(hObject, eventdata, handles); 


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
table_data = get(hObject,'Data'); 
ii = eventdata.Indices; 
% eventdata.EditData
table_data{ii(1),ii(2)} = eventdata.EditData; 
set(hObject, 'Data', table_data);
plotEIs(hObject,handles); 
guidata(hObject,handles); 

function plotEIs(hObject,handles)

table_data = get(handles.uitable1,'Data'); 
axes(handles.axes4); cla; 
hold on; 

positions = handles.datarun.ei.position;
scatter(positions(:,1),positions(:,2),3.5); axis image; axis off; 

if get(handles.showElecs,'Value')
    for e = 1:512
        text(positions(e,1),positions(e,2),num2str(e),'HorizontalAlignment','center'); 
    end
end
colorsneeded = lines(size(table_data,1)); 
for n = 1:1:size(table_data,1)
    if table_data{n,5}
        cellID = table_data{n,1}; 
        cellIndex = get_cell_indices(handles.datarun, cellID);
        ei = handles.datarun.ei.eis{cellIndex}'; % squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        scatter(positions(:,1),positions(:,2),3*eiAmps+0.1,colorsneeded(n,:),'filled');
        elec = table_data{n,2}; 
        text(positions(elec,1),positions(elec,2),num2str(cellID)); 
    end
end


function datapath_Callback(hObject, eventdata, handles)
% hObject    handle to datapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datapath as text
%        str2double(get(hObject,'String')) returns contents of datapath as a double


% --- Executes during object creation, after setting all properties.
function datapath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirname = uigetdir('/Volumes/Analysis','Select Directory'); 
dirname =[dirname dirname(find(dirname == filesep,1,'last'):end)];
set(handles.datapath,'String',dirname)
update_Callback(hObject, eventdata, handles); 

% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ii = eventdata.Indices; 
if isempty(ii)
    return; 
end

% Plot firing rate
axes(handles.axes2);
plot(handles.binranges,handles.firingRates(ii(1),:)); 
ylabel('Spike Rate (Hz)'); xlabel('Time (s)');

%Plot sta
axes(handles.axes1);
sta = handles.datarun.stas.stas{ii(1)}; 
m = max(abs(sta(:)));
a = find(m == abs(sta),1); 
[~, ~, ~, i4] = ind2sub(size(sta), a); 
imagedata = squeeze(sta(:,:,:,i4));
cmax = 0.2; 
cmin = -0.2;  
if size(sta,3) > 1
    imagedata(find(imagedata>=cmax)) = cmax-0.05; 
    imagedata(find(imagedata<=cmin)) = cmin+0.05; 
    imagedata(:,:,1) = (imagedata(:,:,1)-cmin)/(cmax-cmin);
    imagedata(:,:,2) = (imagedata(:,:,2)-cmin)/(cmax-cmin);
    imagedata(:,:,3) = (imagedata(:,:,3)-cmin)/(cmax-cmin);
end
imagesc(imagedata); axis off; colormap gray;

% Plot ei
axes(handles.axes3);
table_data = get(handles.uitable1,'Data'); 
cellID = table_data{ii(1),1}; 
cellIndex = get_cell_indices(handles.datarun, cellID);
ei = handles.datarun.ei.eis{cellIndex}'; % squeeze(ei(1,2:end,:))';
eiAmps = max(ei)-min(ei);
positions = handles.datarun.ei.position;
scatter(positions(:,1),positions(:,2),3*eiAmps+0.1,'filled'); 
axis image;
axis off; 
% data= get(hObject,'Data');
% cell_id = data(ii(1),1); 



% --- Executes on key press with focus on uitable1 and none of its controls.
function uitable1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
disp(eventdata.Character)


% --- Executes on button press in showElecs.
function showElecs_Callback(hObject, eventdata, handles)
% hObject    handle to showElecs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showElecs
plotEIs(hObject,handles); 


% --- Executes on button press in select_all.
function select_all_Callback(hObject, eventdata, handles)
% hObject    handle to select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of select_all


% --- Executes on button press in contourmap.
function contourmap_Callback(hObject, eventdata, handles)
% hObject    handle to contourmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eis = handles.datarun.ei.eis; 
sum = zeros(size(eis{1},1),1); 
table_data = get(handles.uitable1,'Data');
for n = 1:1:size(table_data,1)
    cellID = table_data{n,1};
    cellIndex = get_cell_indices(handles.datarun,cellID); 
    ei = handles.datarun.ei.eis{cellIndex}';
    eiAmps = max(ei) - min(ei); 
    sum = sum + eiAmps';
end
figure; axis image; surf(ei2matrix(sum)); 
figure; axis image; contourf(ei2matrix(sum),24); axis ij; 


% --- Executes on button press in deselectall.
function deselectall_Callback(hObject, eventdata, handles)
% hObject    handle to deselectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    table_data = get(handles.uitable1,'Data');
    
    table_data{:,5} = 0;
    set(handles.uitable1, 'Data', table_data);
end


% --- Executes on button press in axontraces.
function axontraces_Callback(hObject, eventdata, handles)
% hObject    handle to axontraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = figure; 
[xc, yc] = getElectrodeCoords512();
table_data = get(handles.uitable1,'Data');
colors = lines(size(table_data,1)); 
for n = 1:1:size(table_data,1)
    cellID = table_data{n,1};
    cellIndex = get_cell_indices(handles.datarun, cellID);
    ei = handles.datarun.ei.eis{cellIndex}'; % squeeze(ei(1,2:end,:))';
    eiAmps = max(ei)-min(ei);
    thresh = 4;
    [~,col,~] = find(eiAmps > thresh);
    aa = round(eiAmps(col))';
    yy = xc(col)';
    xx = yc(col)';
    
    wy = []; wx = [];
    % dumb way
    for i = 1:1:length(aa)
        tmp = repmat(yy(i),aa(i),1);
        wy = [wy; tmp];
        tmp = repmat(xx(i),aa(i),1);
        wx = [wx; tmp];
    end

    % vector of 1-D look-up table "x" points
    XI = linspace(min(xx),max(xx),max(round(abs(min(xx)-max(xx))/120),2));
    
    % obtain vector of 1-D look-up table "y" points
    YI = lsq_lut_piecewise( wx, wy, XI );
    sortaa = sort(aa,1,'descend'); 
    largestAmps = sortaa(1:2); 
    [~,IA,~] = intersect(aa,largestAmps); 
    
    COMx = 1/sum(largestAmps) * sum(xx(IA).*aa(IA)); 
    COMy = 1/sum(largestAmps) * sum(yy(IA).*aa(IA)); 
       
%     figure(fh); plot(XI,YI,'*-','Color',colors(n,:)); 
% %     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis
% %     hold on; scatter(xx(IA),yy(IA),aa(IA)*6,colors(n,:),'filled'); % largest signals
%     hold on; scatter(COMx,COMy,6*mean(largestAmps), colors(n,:),'filled');
%     hold on; plot(XI,YI,'*-'); 
    
    figure(fh); plot(YI,XI,'*-','Color',colors(n,:)); 
%     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis
%     hold on; scatter(xx(IA),yy(IA),aa(IA)*6,colors(n,:),'filled'); % largest signals
    hold on; scatter(COMy,COMx,6*mean(largestAmps), colors(n,:),'filled');
    text(double(COMy),double(COMx),num2str(cellID)); 
    
%     figure; 
%     hold on; scatter(xc(col),yc(col),eiAmps(col)*6,colors(n,:),'filled');   % Plot eis
%     plot(YI,XI,'*-','Color',0.5*colors(n,:)); 
%     hold on; scatter(COMy,COMx,6*mean(largestAmps), 0.5 * colors(n,:),'filled');
%     text(double(COMy),double(COMx),num2str(cellID)); 
end
axis image; axis off; 
hold on; scatter(xc,yc,5,'black','filled'); 
% hold on; scatter(xc(elecs),yc(elecs),50,'red','filled'); 