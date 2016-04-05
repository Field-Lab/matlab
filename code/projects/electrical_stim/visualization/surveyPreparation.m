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
 
% Last Modified by GUIDE v2.5 24-Aug-2015 12:36:35
 
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
%   Indices: row and column indices of the cell(s) edited
%   PreviousData: previous data for the cell(s) edited
%   EditData: string(s) entered by the user
%   NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%   Error: error string when failed to convert EditData to appropriate value for Data
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
%   Indices: row and column indices of the cell(s) currently selecteds
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
%   Key: name of the key that was pressed, in lower case
%   Character: character interpretation of the key(s) that was pressed
%   Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
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
 
    table_data = get(handles.uitable1,'Data');
   
    
    for r = 1:size(table_data, 1)
        table_data{r,5} = 1;
    end    
    set(handles.uitable1, 'Data', table_data); 
 
 
 
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
[xc,yc] = getElectrodeCoords512();
xq = min(xc):30:max(xc); 
yq = min(yc):30:max(yc); 
vq = griddata(xc,yc,sum,xq,yq');

figure; 
subplot(3,1,1); imagesc(xq,yq,vq); axis xy; axis image;
colorbar; title('interpolated sum of all EIs'); 
subplot(3,1,2); imagesc(xq,yq,log(vq)); axis xy; axis image;
colorbar; title('log(interpolated sum of all EIs)'); 
subplot(3,1,3);  contourf(vq,24); axis xy; axis image; axis off; colorbar; 
title('contour of the sum of EIs'); 
 
 
% --- Executes on button press in deselectall.
function deselectall_Callback(hObject, eventdata, handles)
% hObject    handle to deselectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 
    table_data = get(handles.uitable1,'Data');
   
    
    for r = 1:size(table_data,1)
        table_data{r,5} = 0;
    end    
    set(handles.uitable1, 'Data', table_data);
 
 
 
% --- Executes on button press in lin_axon_traces.
function poly_axon_traces_Callback(hObject, eventdata, handles)
% hObject    handle to lin_axon_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = figure; 
% [xc, yc] = getElectrodeCoords512();
positions = handles.datarun.ei.position;
xc = positions(:,1)'; 
yc = positions(:,2)'; 
table_data = get(handles.uitable1,'Data');
colors = lines(size(table_data,1));
nearby_axons = zeros(1, size(xc,2));
nearby_somas = zeros(1, size(xc,2));
nearby_range = 1; %measured in number of electrode distances, can be fractional
 
% load('validIDs.mat');
 
for n = 1:1:size(table_data,1)
    cellID = table_data{n,1};
    cellIndex = get_cell_indices(handles.datarun, cellID);
    ei = handles.datarun.ei.eis{cellIndex}'; % squeeze(ei(1,2:end,:))';
    eiAmps = max(ei)-min(ei);
    thresh = 5;
    [~,col,~] = find(eiAmps > thresh);
    aa = round(eiAmps(col))';
   
    sortaa = sort(aa,1,'descend'); 
    largestAmps = sortaa(1:2); 
  
    [XI, YI, ~, COMx, COMy, valid] = weighted_axon_poly_reg(eiAmps);
    
    figure(fh); 
    if valid || cellID == 2796 || cellID == 2842 || cellID == 2626 || cellID == 3173
        plot(XI,YI,'-','Color',colors(n,:));
        hold on;   
        scatter(COMx,COMy,mean(largestAmps), colors(n,:),'filled');
    else
        fprintf('(above warning for axon %d)\n',cellID);  
    end 
    
    
%    text(double(COMx),double(COMy),num2str(cellID)); 
%     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis 
    
    close = zeros(size(xc,2), 1);
    
    if valid || cellID == 2796 || cellID == 2842 || cellID == 2626 || cellID == 3173
    %Cuts down on axon steps to reduce runtime
    XI = XI(1:floor(length(XI)/10):length(XI));
    YI = YI(1:floor(length(YI)/10):length(YI));
    
    %Finding axons within one electrode spacing
    for elec = 1:size(xc,2)
        for ind = 1:size(YI,2)
            if pdist([XI(ind) YI(ind); xc(elec), yc(elec)]) < (60 * nearby_range)
                close(find(close == 0, 1, 'first')) = elec;
                break;
            end
        end
    end
    close(close == 0) = [];
    close = unique(close);
    nearby_axons(close) = nearby_axons(close) + 1;
    end
    
    %Finding somas within one elctrode spacing
    close = zeros(size(xc,2), 1);
    for ind = 1:size(xc, 2)
        if pdist([COMx COMy; xc(ind) yc(ind)]) < (60 * nearby_range)
            close(find(close == 0, 1, 'first')) = ind;
        end
    end
    close(close == 0) = [];
    close = unique(close);
    nearby_somas(close) = nearby_somas(close) + 1;
 
%     hold on; scatter(xx(IA),yy(IA),aa(IA)*6,colors(n,:),'filled'); % largest signals
   
    
%     figure; 
%     hold on; scatter(xc(col),yc(col),eiAmps(col)*6,colors(n,:),'filled');   % Plot eis
%     plot(YI,XI,'*-','Color',0.5*colors(n,:)); 
%     hold on; scatter(COMy,COMx,6*mean(largestAmps), 0.5 * colors(n,:),'filled');
%     text(double(COMy),double(COMx),num2str(cellID)); 
end
axis image; axis off; 
hold on; scatter(xc,yc,5,'black','filled'); 
 
figure; scatter(xc,yc,300,nearby_axons,'filled'); colorbar; title('Axons');
xlabel(['# axons within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
figure; scatter(xc,yc,300,nearby_somas,'filled'); colorbar; title('Somas');
xlabel(['# somas within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
figure; 
for x = 1:512
    text(xc(x)+20,yc(x)+20,num2str(x),'HorizontalAlignment','center', 'Color', 'white');
end
hold on;
scatter(xc,yc,300,nearby_somas+nearby_axons,'filled'); colorbar; title('Both');
xlabel(['# axons+somas within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
handles.poly_nearby_axons = nearby_axons;
handles.poly_nearby_somas = nearby_somas;
guidata(hObject, handles);
 
 
 
% --- Executes on button press in linGenNearbyFeaturesList.
function linGenNearbyFeaturesList_Callback(hObject, eventdata, handles)
% hObject    handle to linGenNearbyFeaturesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Axons = handles.lin_nearby_axons';
Somas = handles.lin_nearby_somas';
Both = Axons+Somas;
 
sortby = Both;
 
[sorted, Electrodes] = sort(sortby, 1, 'descend');
 
Axons = Axons(Electrodes);
Somas = Somas(Electrodes);
Both = Axons+Somas;
 
sorted_full = [Electrodes Axons Somas Both];
cnames = {'Electrode', 'Axons', 'Somas', 'Both'};
t = uitable(figure(), 'Data', sorted_full, 'ColumnName', cnames);
 
t.Position(3) = t.Extent(3);
%table(Electrodes, Axons, Somas, Both);
 
 
% --- Executes on button press in poly_axon_traces.
function lin_axon_traces_Callback(hObject, eventdata, handles)
fh = figure; 
% [xc, yc] = getElectrodeCoords512();
positions = handles.datarun.ei.position;
xc = positions(:,1)'; 
yc = positions(:,2)'; 
table_data = get(handles.uitable1,'Data');
colors = lines(size(table_data,1));
nearby_axons = zeros(size(xc));
nearby_somas = zeros(size(xc));
nearby_range = 1; %measured in number of electrode distances, can be fractional
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
    close = zeros(512, 1);
    for ind = 1:(size(YI, 1)-1)
        p1 = [YI(ind) XI(ind)]; p2 = [YI(ind+1) XI(ind+1)];
        vec = [-(p1(2)-p2(2)) p1(1)-p2(1)];
        vec = vec/norm(vec)*nearby_range*60;
        c1 = p1+vec;
        c2 = p1-vec;
        c3 = p2-vec;
        c4 = p2+vec;
        xverts = [c1(1) c2(1) c3(1) c4(1)]; yverts = [c1(2) c2(2) c3(2) c4(2)];
        contained = find(inpolygon(xc, yc, xverts, yverts));
        
        for elec = contained
            if l2pdist([p1; p2],[xc(elec), yc(elec)]) < 60
                close(find(close == 0, 1, 'first')) = elec;
            end
        end
    end
    close(close == 0) = [];
    close = unique(close);
    nearby_axons(close) = nearby_axons(close) + 1;
    close = zeros(512, 1);
    for ind = 1:size(xc, 2)
        if pdist([COMy COMx; xc(ind) yc(ind)]) < 60
            close(find(close == 0, 1, 'first')) = ind;
        end
    end
    close(close == 0) = [];
    close = unique(close);
    nearby_somas(close) = nearby_somas(close) + 1;
    
 
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
%hold on; scatter(xc,yc,5,'black','filled'); 
 
figure; scatter(xc,yc,300,nearby_axons,'filled'); colorbar; title('Axons');
xlabel(['# axons within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
figure; scatter(xc,yc,300,nearby_somas,'filled'); colorbar; title('Somas');
xlabel(['# somas within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
figure; 
for x = 1:512
    text(xc(x)+20,yc(x)+20,num2str(x),'HorizontalAlignment','center', 'Color', 'white');
end
hold on;
scatter(xc,yc,300,nearby_somas+nearby_axons,'filled'); colorbar; title('Both');
xlabel(['# axons+somas within ' num2str(nearby_range) ' elec distance(s)']);
axis image; axis off;
set(findall(gca, 'type', 'text'), 'visible', 'on');
handles.lin_nearby_axons = nearby_axons;
handles.lin_nearby_somas = nearby_somas;
guidata(hObject, handles);
 
% hold on; scatter(xc(elecs),yc(elecs),50,'red','filled'); 
 
 
% --- Executes on button press in polyGenNearbyFeaturesList.
function polyGenNearbyFeaturesList_Callback(hObject, eventdata, handles)
% hObject    handle to linGenNearbyFeaturesList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Axons = handles.poly_nearby_axons';
Somas = handles.poly_nearby_somas';
Both = Axons+Somas;
 
sortby = Both;
 
[sorted, Electrodes] = sort(sortby, 1, 'descend');
 
Axons = Axons(Electrodes);
Somas = Somas(Electrodes);
Both = Axons+Somas;
 
sorted_full = [Electrodes Axons Somas Both];
cnames = {'Electrode', 'Axons', 'Somas', 'Both'};
t = uitable(figure(), 'Data', sorted_full, 'ColumnName', cnames);
 
t.Position(3) = t.Extent(3);
 
 
% --- Executes on button press in stim_bundles.
function stim_bundles_Callback(hObject, eventdata, handles)
% hObject    handle to stim_bundles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ei_path = get(handles.datapath,'String');
slash = strfind(ei_path,filesep);
main_dir = ei_path(1:max(slash)-1);
 
input = inputdlg('Enter number of movie folders');
num_folders = str2num(input{1});
folders = cell(1,num_folders);
 
for f=1:num_folders;
    dataPath = uigetdir(main_dir,strcat('Select directory #',num2str(f)));
    folders{1,f} = dataPath;
end
 
% [elec_x, elec_y] = getElectrodeCoords512();
positions = handles.datarun.ei.position;
elec_x = positions(:,1)'; 
elec_y = positions(:,2)'; 
end_amps = zeros(512,1);
 
for n = 1:num_folders
    files = dir(folders{n});
    files = {files(cell2mat({files(:).isdir})).name};

    for d = 1:length(files)
        if length(files{d}) > 1 && strcmp(files{d}(1),'p') && isstrprop(files{d}(2), 'digit')
            elec_num = str2num(files{d}(2:length(files{d})));
            movieNos = findMovieNos(folders{n}, elec_num);
            [rawData, ~] = generateEiFromStimPattern(folders{n}, elec_num,'movieNo',max(movieNos), 'suppressPlots', true);
             avg_trials = squeeze(mean(rawData,1));
             %10 to 40 are the timepoints that avoid the artifact. May need
             %adjustment for different experiments
             after_artifact = avg_trials(:,10:40);
             spikes = max(after_artifact, [], 2) - min(after_artifact, [], 2);
             end_amps = end_amps + spikes;
        end    
    end
     
end
    figure;
    scatter(elec_x,elec_y,350,end_amps/512,'filled');
    axis off; axis image; c = colorbar; caxis([0 max(end_amps/512)]);
     ylabel(c,'  \muV','rot',0); set(gca,'FontSize',16);

