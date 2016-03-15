function varargout = clusterSorting(varargin)
% CLUSTERSORTING MATLAB code for clusterSorting.fig
%      CLUSTERSORTING, by itself, creates a new CLUSTERSORTING or raises the existing
%      singleton*.
%
%      H = CLUSTERSORTING returns the handle to a new CLUSTERSORTING or the handle to
%      the existing singleton*.
%
%      CLUSTERSORTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTERSORTING.M with the given input arguments.
%
%      CLUSTERSORTING('Property','Value',...) creates a new CLUSTERSORTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before clusterSorting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to clusterSorting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help clusterSorting

% Last Modified by GUIDE v2.5 11-Mar-2016 16:25:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clusterSorting_OpeningFcn, ...
                   'gui_OutputFcn',  @clusterSorting_OutputFcn, ...
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


% --- Executes just before clusterSorting is made visible.
function clusterSorting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to clusterSorting (see VARARGIN)

% Choose default command line output for clusterSorting
handles.output = hObject;

% Populate Electrode Array
axis(handles.axes0);
% handles.axes0 = axes('NextPlot','replacechildren');
[handles.xCoords, handles.yCoords] = getElectrodeCoords512(); 
handles.array = scatter(handles.xCoords,handles.yCoords,500,'filled','black');
axis image; axis off; 
for i = 1:length(handles.xCoords)
    text(handles.xCoords(i),handles.yCoords(i)-30,num2str(i),'HorizontalAlignment','center','Color','black'); 
end
set(handles.array,'ButtonDownFcn',{@arrayClickCallback, handles}); 
caxis([0 4]); %default caxis
% handles.cbar = colorbar('EastOutside','YTickLabel',' '); 
handles.cbar = colorbar('EastOutside');
handles.map = flipud(hot); colormap(handles.map); 
set(hObject,'CloseRequestFcn',@my_closereq)
set(handles.dataPanel,'Visible','off'); 
set(handles.slider1,'SliderStep',[1 1]); 
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes clusterSorting wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function arrayClickCallback (hObject, eventdata, localhandles)

axesHandle = get(hObject,'Parent');
latestHandles = guidata(axesHandle); 
a = get(axesHandle,'CurrentPoint');
x = a(1,1);
y = a(1,2);
diffX = localhandles.xCoords - x;
diffY = localhandles.yCoords - y;
d = hypot(diffX,diffY); 
pattern = find(d == min(d));
set(localhandles.patternNo,'String',num2str(pattern)); 
set(localhandles.templateElectrode,'String',num2str(pattern)); 
% Check activationResults.mat
dirname = get(latestHandles.datapath,'String');
resultsFile = fullfile(dirname,'activationResults.mat');
if ~exist(resultsFile,'file')
    pattern = struct('thresholds',4*ones(512,1),'activeNeuronID',[],'userAccepted',zeros(512,1)); 
    activationResults = struct('pattern',repmat(pattern,512,1)); 
    save(resultsFile,'activationResults');
    
else
    tmp = load(resultsFile); 
    activationResults = tmp.activationResults; 
end
latestHandles = go_Callback(hObject,eventdata,latestHandles);
latestHandles.test = 'arrayClickCallback can update handles structure'; 
guidata(axesHandle,latestHandles); % Update latest handles structure (not local copy)
disp(['pattern number: ' num2str(pattern)]);

% --- Outputs from this function are returned to the command line.
function varargout = clusterSorting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
set(handles.datapath,'String',dirname)
% Check for activation results file
resultsFile = fullfile(dirname,'activationResults.mat');
if ~exist(resultsFile,'file')
    pattern = struct('thresholds',4*ones(512,1),'activeNeuronID',[],'userAccepted',zeros(512,1)); 
    activationResults = struct('pattern',repmat(pattern,512,1)); 
    save(resultsFile,'activationResults');
end
guidata(hObject,handles); 
% update_Callback(hObject, eventdata, handles); 


function patternNo_Callback(hObject, eventdata, handles)
% hObject    handle to patternNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of patternNo as text
%        str2double(get(hObject,'String')) returns contents of patternNo as a double
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end
handles = go_Callback(hObject,eventdata,handles); 
guidata(hObject,handles); 

% --- Executes during object creation, after setting all properties.
function patternNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patternNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes0_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('test matlab auto version');

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    % Increase pattern number
    patternNo = str2double(get(handles.patternNo,'String')); 
    if patternNo < 512
        patternNo = patternNo + 1; 
    end
    set(handles.patternNo,'String',num2str(patternNo)); 
    patternNo_Callback(hObject,eventdata,handles); 
    set(hObject,'Value',0.5);
else
    % Decrease pattern number
    patternNo = str2double(get(handles.patternNo,'String')); 
    if patternNo > 1
        patternNo = patternNo - 1;
    end
    set(handles.patternNo,'String',num2str(patternNo)); 
    patternNo_Callback(hObject,eventdata,handles); 
    set(hObject,'Value',0.5);
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in go.
function handles = go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for ii = 1:6
    command = sprintf(['cla(handles.axes%d,''reset'');'...
        'axis(handles.axes%d,''off'');'...
        'set(handles.ignore%d,''Value'',0);'],ii,ii,ii);
    eval(command);
end
if isfield(handles,'tempData')
    handles = rmfield(handles,'tempData');
end
if isfield(handles,'vectorMultiplier')
    handles = rmfield(handles,'vectorMultiplier'); 
end
if isfield(handles,'responseRate')
    handles = rmfield(handles,'responseRate'); 
end
guidata(hObject,handles); 
pathToAnalysisData = get(handles.datapath,'String'); 
patternNo = str2double(get(handles.patternNo,'String'));
sampleRange = 5:35;
movieNos = findMovieNos(pathToAnalysisData,patternNo);

expU = 2; %exponent for partition matrix U
maxIt = 8; % maximum number of iterations
minImp = 1e-5; %minimum amount of improvement
fcmVerbose = 0; %Information displayed during iteration
clusteringOptions = [expU maxIt minImp fcmVerbose];
mThresh = 0.85; %membership function threshold
numAllowedOutliers = 2; %number of outliers allowed
minTracesPerCluster = 3; % minimum number of traces per cluster required to register as activation
plotOn = 0;

tempRecElec = patternNo;
fprintf('Calculating results at electrode %d stimulating with electrode %d\n',tempRecElec,patternNo); 
% figure; set(gcf,'Position',[15 145 1600 960]);
set(handles.dataPanel,'Visible','on'); 
subplotIdx = 1; 
for m = 1:size(movieNos,2)
    % Load stim amps
    stimAmp = getStimAmps(pathToAnalysisData, patternNo, movieNos(m), 'numElectrodes', 512);
    
    % Load data
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(m), 99999);
    testTrace = squeeze(dataTraces(:,tempRecElec,sampleRange));
    subMeanTrace = testTrace - repmat(mean(testTrace),size(testTrace,1),1);
    if get(handles.useElectrodeCluster,'Value')
        clusterRecElecs = getCluster512(tempRecElec);
        testTraceCluster = squeeze(dataTraces(:,clusterRecElecs,sampleRange));
        subMeanTraceCluster = testTraceCluster - repmat(mean(testTraceCluster),size(testTraceCluster,1),1);
        clusterTraces = reshape(permute(subMeanTraceCluster,[1 3 2]),size(dataTraces,1),[]);
        subMeanTrace = clusterTraces; % hack to test using surrounding electrodes
    end
    
    % Do PCA
    [~,score] = pca(subMeanTrace);
    
    % Fuzzy C-means clustering on the first 2 principal components
    fcmdata = score(:,1:2);
    [~, U, ~] = fcm(fcmdata, str2double(get(handles.numClusters,'String')), clusteringOptions);
    maxU = max(U);
    index1 = find(U(1,:) == maxU);
    index2 = find(U(2,:) == maxU);

    if length(index2) >= minTracesPerCluster && length(index1) >= minTracesPerCluster
        % Use values of the membership function to determine if clustering
        % was strong, this is the metric for finding bifurcations. 
        if  (sum(U(1,index1) > mThresh) >= (length(index1)-numAllowedOutliers)) || (sum(U(2,index2) > mThresh) >= (length(index2)-numAllowedOutliers))
            set(handles.dataPanel,...
                'Title',sprintf('stimulating electrode %0.0f, recording electrode %0.0f',patternNo,tempRecElec));  
            plotOn = 1;
        end
    end
    
    while plotOn
        if any(sum(abs(subMeanTrace))<20) && get(handles.useElectrodeCluster,'Value')==0 % Saturation not as important if using cluster of elecs
            plotOn = 0;
            try
                surprise = getEmoji();
                msgbox(sprintf('movie %d is saturated for pattern %d at electrode %d',movieNos(m),patternNo,tempRecElec),'SATURATED SIGNAL','custom',surprise);
            catch
                msgbox(sprintf('movie %d is saturated for pattern %d at electrode %d',movieNos(m),patternNo,tempRecElec),'SATURATED SIGNAL');
            end
            handles.responseRate = responseRate; 
            handles.responseStimAmps = responseStimAmps; 
            handles.vectorMultiplier = ones(size(responseRate)); 
            handles.tempRecElec = tempRecElec; 
            guidata(hObject,handles);
%             calculateRespRate_Callback(hObject, eventdata, handles); 
            return;
        end
        
        group1 = testTrace(index1,:);
        group2 = testTrace(index2,:);
        
        eval(['axes(handles.axes' num2str(subplotIdx) ');']);
        
        if min(mean(group1,1)) < min(mean(group2,1))
            plot(subMeanTrace(index2,:)','Color',0.7*[1 1 1]); hold on;
            plot(subMeanTrace(index1,:)','red');
            axis off;
            title(sprintf('movie %d',movieNos(m)));
            handles.tempData(subplotIdx).spikes = subMeanTrace(index1,:)';
            handles.tempData(subplotIdx).misses = subMeanTrace(index2,:)';
            
            responseStimAmps(subplotIdx) = stimAmp;
            responseRate(subplotIdx) = length(index1)/size(testTrace,1);
        else
            plot(subMeanTrace(index1,:)','Color',0.7*[1 1 1]); hold on;
            plot(subMeanTrace(index2,:)','red');
            axis off;
            title(sprintf('movie %d',movieNos(m)));
            handles.tempData(subplotIdx).spikes = subMeanTrace(index2,:)';
            handles.tempData(subplotIdx).misses = subMeanTrace(index1,:)';
            responseStimAmps(subplotIdx) = stimAmp;
            responseRate(subplotIdx) = length(index2)/size(testTrace,1);
        end
        
       plotOn = 0;
        subplotIdx = subplotIdx + 1;
        if subplotIdx>6 
            
            handles.responseRate = responseRate; 
            handles.responseStimAmps = responseStimAmps; 
            handles.vectorMultiplier = ones(size(responseRate)); 
            handles.tempRecElec = tempRecElec; 
            guidata(hObject,handles);
           
            if get(handles.displayCurrentElec,'Value');
                setElecColor(handles)
            elseif get(handles.allElectrodeDisplay,'Value');
                setElecColorShowPast(handles);
            end
            return;
        end
    end
    plotOn = 0;
end
if exist('responseRate','var')
    handles.responseRate = responseRate;
    handles.responseStimAmps = responseStimAmps;
else
    handles.responseRate = 0;
    handles.responseStimAmps = 1;
end
handles.vectorMultiplier = ones(size(handles.responseRate));
handles.tempRecElec = tempRecElec;
guidata(hObject,handles);

function setElecColor(handles)
% Check for activation results file
dirname = get(handles.datapath,'String');
resultsFile = fullfile(dirname,'activationResults.mat');
if exist(resultsFile,'file')
    tmp = load(resultsFile); 
    activationResults = tmp.activationResults; 
end
patternNo = str2double(get(handles.patternNo,'String')); 
info = activationResults.pattern(patternNo); 
toPlot = info.thresholds.*info.userAccepted; 
maxAmp = str2double(get(handles.maxAmp,'String')); 
minAmp = str2double(get(handles.minAmp,'String'));  
axes(handles.axes0); hold on; 
h = scatter(handles.xCoords, handles.yCoords, 500, toPlot,'filled'); 
handles.map(1,:) = 0.5*[1 1 1]; 
colormap(handles.map); 
axis off; axis image;   
caxis([minAmp maxAmp]);
ylabel(handles.cbar,'  \muA','rot',0);
set(gca,'FontSize',16);
set(h,'ButtonDownFcn',{@arrayClickCallback, handles});

function setElecColorShowPast(handles)
% Check for activation results file
dirname = get(handles.datapath,'String');
resultsFile = fullfile(dirname,'activationResults.mat');
if exist(resultsFile,'file')
    tmp = load(resultsFile); 
    activationResults = tmp.activationResults; 
end
toPlot = zeros(512,1); 
for e = 1:512
    toPlot(e) = activationResults.pattern(e).thresholds(e);
end

maxAmp = str2double(get(handles.maxAmp,'String')); 
minAmp = str2double(get(handles.minAmp,'String')); 
axes(handles.axes0); hold on; 
h = scatter(handles.xCoords, handles.yCoords, 500, toPlot,'filled'); 
handles.map(1,:) = 0.5*[1 1 1]; 
colormap(handles.map); 
axis off; axis image;   
caxis([minAmp maxAmp]);
ylabel(handles.cbar,'  \muA','rot',0);
set(gca,'FontSize',16);
set(h,'ButtonDownFcn',{@arrayClickCallback, handles});

function showTemplates(tempRecElec, sampleRange, templatepath)

% Project onto template waveform
dirname = templatepath; 
neuronPath =[dirname dirname(find(dirname == filesep,1,'last'):end) '.neurons'];
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
neuronIds = neuronFile.getIDList(); %returns a list of all the neuron IDs in this neuron file


ei_path =[dirname dirname(find(dirname == filesep,1,'last'):end) '.ei'];
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(ei_path);

neuronId = neuronIds(1);
neuronEI = eiFile.getImage(neuronId);
neuronEI_volt = squeeze(neuronEI(1,:,:)).'; % Actual EI
eiMatrix = zeros(length(neuronIds),size(neuronEI_volt,1),512);

% LOAD ONE EI AT A TIME
% loop through each cell
for cc = 1:length(neuronIds)
    neuronId = neuronIds(cc);
    neuronEI = eiFile.getImage(neuronId);
    neuronEI_volt = squeeze(neuronEI(1,:,:)).'; % Actual EI
    
    % throw out the triggers
    ei = neuronEI_volt(:,2:end);
    % restore double
    allEis{neuronId} = double(ei);
    eiMatrix(cc,:,:) = double(ei);
end

% Check waveforms on electrode of interest.
waveformsOn1Elec = squeeze(eiMatrix(:,:,tempRecElec));
[rowIdx,~]= find(waveformsOn1Elec<-30);
% templates = waveformsOn1Elec(unique(rowIdx),sampleRange);
clusterElecs = getCluster512(tempRecElec); 
figure; set(gcf, 'Position', [60 956 1670 150]); 
for i = 1:length(clusterElecs)
    subplot(1,length(clusterElecs),i);
    waveformsOn1Elec = squeeze(eiMatrix(:,:,clusterElecs(i)));
    templates = waveformsOn1Elec(unique(rowIdx),sampleRange);
    plot(templates','LineWidth',2);
    title(num2str(clusterElecs(i)));
end
legend(num2str(neuronIds(unique(rowIdx))));


% --- Executes on button press in templateDisp.
function templateDisp_Callback(hObject, eventdata, handles)
% hObject    handle to templateDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tempRecElec = str2double(get(handles.templateElectrode,'String'));
if isnan(tempRecElec)
    tempRecElec = str2double(get(handles.patternNo,'String'));
end
sampleRange        = 1:40;
templatepath = get(handles.templatepath,'String');
showTemplates(tempRecElec, sampleRange, templatepath);


function templateElectrode_Callback(hObject, eventdata, handles)
% hObject    handle to templateElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of templateElectrode as text
%        str2double(get(hObject,'String')) returns contents of templateElectrode as a double


% --- Executes during object creation, after setting all properties.
function templateElectrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to templateElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function my_closereq(hObject, eventdata, handles)  
% User-defined close request function
% to display a question dialog box

selection = questdlg('Close Cluster Sorting GUI?',...
                     'Close Request Function',...
                     'Yes','Nooo','Yes');
switch selection,
   case 'Yes',
     delete(hObject);
   case 'Nooo'
     return
end


% --- Executes on button press in ignore1.
function ignore_Callback(hObject, eventdata, handles)
% hObject    handle to ignore1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
latestHandle = guidata(get(hObject,'Parent')); 
tag =  get(hObject,'Tag');
idx = str2double(tag(end)); 
if get(hObject,'Value')
    handles.vectorMultiplier(idx) = 0;
else
    handles.vectorMultiplier(idx) = 1;
end
guidata(hObject,handles); 

% --- Executes on button press in calculateRespRate.
function calculateRespRate_Callback(hObject, eventdata, handles)
% hObject    handle to calculateRespRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = guidata(get(hObject,'Parent')); 
responseStimAmps = nonzeros(handles.responseStimAmps.*handles.vectorMultiplier)';
responseRate = nonzeros(handles.responseRate.*handles.vectorMultiplier)';
disp(['using the following movies to calculate threshold: '...
    num2str(handles.vectorMultiplier)]); 
[abs([responseStimAmps(1)/1.1 responseStimAmps responseStimAmps(end)*1.1]) ; 0 responseRate 1]
forcedStimAmps = abs([responseStimAmps(1)/1.1 responseStimAmps responseStimAmps(end)*1.1]); 
forcedRespRate = [0 responseRate 1];
% Define function that will be used to fit data
% (F is a vector of fitting parameters)
f = @(F,x) (1 +exp(-F(1)*(x - F(2)))).^(-1); % sigmoid
F_fitted = nlinfit(forcedStimAmps,forcedRespRate,f,[1 1]);

% Plot data fit
stimAmps = forcedStimAmps;
xx = stimAmps(1):0.001:stimAmps(end);
yy = f(F_fitted,xx);

% y = f(F_fitted,abs(responseStimAmps));
projectionComplete(1,:) = xx;
projectionComplete(2,:) = yy;

% [~, projectionComplete, error] = gaussCdfFitter([abs(responseStimAmps); responseRate]);
figure; plot(abs(responseStimAmps), responseRate,'-x');
hold on; plot(projectionComplete(1,:),projectionComplete(2,:),'-g');
idx = find(projectionComplete(2,:) > 0.5, 1, 'first');
actThresh = projectionComplete(1,idx);
hold on; line([projectionComplete(1,1) actThresh],[0.5 0.5],'Color','r');
if actThresh
    hold on; line([actThresh actThresh],[projectionComplete(2,1) 0.5],'Color','r');
else
    actThresh = 4;
end
legend('data','gauss cdf fit');
xlabel('stimulation amplitude (uA)');
ylabel('response rate');
title(sprintf('Electrode no. %d; Approx activation threshold %0.2f uA',handles.tempRecElec,actThresh));
handles.tempActThresh = actThresh; 
guidata(hObject,handles); 

% --- Executes on button press in acceptThreshold.
function acceptThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to acceptThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for activation results file
dirname = get(handles.datapath,'String');
resultsFile = fullfile(dirname,'activationResults.mat');
if ~exist(resultsFile,'file')
    pattern = struct('thresholds',4*ones(512,1),'activeNeuronID',[],'userAccepted',zeros(512,1)); 
    activationResults = struct('pattern',repmat(pattern,512,1)); 
    save(resultsFile,'activationResults'); 
else
    tmp = load(resultsFile); 
    activationResults = tmp.activationResults; 
end
activationThresholds = activationResults.pattern(str2double(get(handles.patternNo,'String'))).thresholds;
activationThresholds(handles.tempRecElec)= handles.tempActThresh;
activationResults.pattern(str2double(get(handles.patternNo,'String'))).thresholds = activationThresholds; 
activationResults.pattern(str2double(get(handles.patternNo,'String'))).userAccepted(handles.tempRecElec) = 1;  %#ok<STRNU>

save(resultsFile,'activationResults'); 
handles.actThresh = handles.tempActThresh; 
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end

function maxAmp_Callback(hObject, eventdata, handles)
% hObject    handle to maxAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','4'); 
end
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end 

% --- Executes during object creation, after setting all properties.
function maxAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minAmp_Callback(hObject, eventdata, handles)
% hObject    handle to minAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','0'); 
end
% callback to refresh electrode map.
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end

% --- Executes during object creation, after setting all properties.
function minAmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in switchgroups1.
function switchgroups_Callback(hObject, eventdata, handles)
% hObject    handle to switchgroups1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tag =  get(hObject,'Tag');
idx = str2double(tag(end));
if idx <= size(handles.tempData,2)
    spikes = handles.tempData(idx).misses;
    misses = handles.tempData(idx).spikes;
    
    eval(sprintf('cla(handles.axes%d); axes(handles.axes%d);',idx,idx));
    plot(misses,'Color',0.7*[1 1 1]); hold on; plot(spikes,'red');
    axis off;
    
    handles.tempData(idx).spikes = spikes;
    handles.tempData(idx).misses = misses;
    
    handles.responseRate(idx) = size(spikes,2)/(size(spikes,2)+size(misses,2));
else
    disp('no distinct spikes / misses');
end
set(hObject,'Value',0);
guidata(hObject,handles);


% --- Executes on button press in subtract1.
function subtract_Callback(hObject, eventdata, handles)
% hObject    handle to subtract1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% axesHandle = get(hObject,'Parent');
% latestHandles = guidata(handles.axes0); 
tag =  get(hObject,'Tag');
idx = str2double(tag(end));
if idx <= size(handles.tempData,2)
    spikes = handles.tempData(idx).spikes;
    misses = handles.tempData(idx).misses;
    figure; plot(mean(spikes,2) - mean(misses,2));
    title('Spikes - misses');
else
    disp('no distinct spikes / misses'); 
end
set(hObject,'Value',0);


% --- Executes on button press in useElectrodeCluster.
function useElectrodeCluster_Callback(hObject, eventdata, handles)
% hObject    handle to useElectrodeCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in manualentry.
function manualentry_Callback(hObject, eventdata, handles)
% hObject    handle to manualentry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = inputdlg('Enter approximate activation threshold',...
    'Manual entry',[1 50]);
manualThreshold = str2double(answer);

% Check for activation results file
dirname = get(handles.datapath,'String');
resultsFile = fullfile(dirname,'activationResults.mat');
if ~exist(resultsFile,'file')
    pattern = struct('thresholds',4*ones(512,1),'activeNeuronID',[],'userAccepted',zeros(512,1)); 
    activationResults = struct('pattern',repmat(pattern,512,1)); 
    save(resultsFile,'activationResults'); 
else
    tmp = load(resultsFile); 
    activationResults = tmp.activationResults; 
end
activationThresholds = activationResults.pattern(str2double(get(handles.patternNo,'String'))).thresholds;
activationThresholds(handles.tempRecElec)= manualThreshold;
activationResults.pattern(str2double(get(handles.patternNo,'String'))).thresholds = activationThresholds; 
activationResults.pattern(str2double(get(handles.patternNo,'String'))).userAccepted(handles.tempRecElec) = 1;  %#ok<STRNU>

save(resultsFile,'activationResults'); 
handles.actThresh = handles.tempActThresh; 
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end


function templatepath_Callback(hObject, eventdata, handles)
% hObject    handle to templatepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of templatepath as text
%        str2double(get(hObject,'String')) returns contents of templatepath as a double


% --- Executes during object creation, after setting all properties.
function templatepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to templatepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseEItemplate.
function browseEItemplate_Callback(hObject, eventdata, handles)
% hObject    handle to browseEItemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startpath = get(handles.datapath,'String'); 
% Project onto template waveform
templatePath = uigetdir(startpath,'Select ei directory');
set(handles.templatepath,'String',templatePath); 


% --- Executes on button press in allElectrodeDisplay.
function allElectrodeDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to allElectrodeDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allElectrodeDisplay
if get(hObject,'Value')
    set(handles.displayCurrentElec,'Value',0); 
else
    set(handles.displayCurrentElec,'Value',1); 
end
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end

% --- Executes on button press in displayCurrentElec.
function displayCurrentElec_Callback(hObject, eventdata, handles)
% hObject    handle to displayCurrentElec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayCurrentElec
if get(hObject,'Value')
    set(handles.allElectrodeDisplay,'Value',0); 
else
    set(handles.allElectrodeDisplay,'Value',1); 
end
if get(handles.displayCurrentElec,'Value');
    setElecColor(handles)
elseif get(handles.allElectrodeDisplay,'Value');
    setElecColorShowPast(handles);
end



function numClusters_Callback(hObject, eventdata, handles)
% hObject    handle to numClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numClusters as text
%        str2double(get(hObject,'String')) returns contents of numClusters as a double


% --- Executes during object creation, after setting all properties.
function numClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Enforce 2-4 clusters. 
if str2double(get(hObject,'String'))>4
    set(hObject,'String','2'); 
end
if str2double(get(hObject,'String'))<2
    set(hObject,'String','2'); 
end

% go_Callback(hObject, eventdata, handles)