clear;

rootDataFolder = '/media/MEA_PROCESSED_5/2013-05-23-0/data';
rootLogfilesFolder = '/media/MEA_PROCESSED_5/2013-05-23-0/logfiles/';

% Data files that you want to process
dataStrs = {'002'};
 
for kk=1:length(dataStrs)
    statsFolder = fullfile(rootDataFolder,sprintf('data%s/statistics/activation', dataStrs{kk}));
    figuresFolder = fullfile(rootDataFolder,sprintf('data%s/figures/activation', dataStrs{kk}));
    thresholdsFolder = fullfile(rootDataFolder,sprintf('data%s/statistics/thresholds', dataStrs{kk}));
        
    plotActivation(statsFolder,figuresFolder,thresholdsFolder,...
        'imageFormat','epsc');
end