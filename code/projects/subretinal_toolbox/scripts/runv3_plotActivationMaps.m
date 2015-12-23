rootDataFolder = '/media/MEA_PROCESSED_5/2013-06-20-0/data/';
rootLogfilesFolder = '/media/MEA_PROCESSED_5/2013-06-20-0/logfiles/';

% Setting processing parameters
stimulatedNeurons = [];  

powerToIrradianceFactor = 11.2; % Scales the logfile powers to irradiances

% Data files that you want to process
dataStrs = {'004'};
 
for kk=1:length(dataStrs)
    dataFolder = [rootDataFolder 'data' dataStrs{kk} '/vision_processing/data' dataStrs{kk}];
    statsFolder = [rootDataFolder 'data' dataStrs{kk} '/statistics/activation'];
    figuresFolder = [rootDataFolder 'data' dataStrs{kk} '/figures/maps'];
    logFilePath = [rootLogfilesFolder 'logfile' dataStrs{kk} '.txt'];
    
    [M, paramNames] = readLogFile(logFilePath);
    pulseDuration = M(:,ismember(paramNames,'Pulse Duration'));
    pulsePower = M(:,ismember(paramNames,'Power'));
    
    for ll=1:length(pulsePower)
        titleString = sprintf('Activation map, %dms pulse, %4.2fmW/mm^2 irradiance',...
            pulseDuration(ll),pulsePower(ll)*powerToIrradianceFactor);
        
        plotPositionEI(dataFolder,statsFolder,figuresFolder,...
            'stimulatedNeuronsList',stimulatedNeurons,... 
            'titleString',titleString,...
            'stimToPlot',ll,...       %             'positionPixel',pixelPositionList(ll,:),...
            'savePlot',true,...
            'displayNeuronNames',true,...
            'maxNSpikes',2,...
            'nColorsColormap',10,...
            'imageFormat','epsc');
    end
end
