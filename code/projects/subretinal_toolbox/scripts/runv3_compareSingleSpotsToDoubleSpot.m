clear;

%% 2012-08-06

singleSpotStatsFolder1 = '/media/MEA_PROCESSED_4/2012-08-06-0/data/data001/statistics/activation';
singleSpotStatsFolder2 = singleSpotStatsFolder1;
doubleSpotStatsFolder = '/media/MEA_PROCESSED_4/2012-08-06-0/data/data003/statistics/activation';

doubleSpotRun = [2:8];
singleSpotRun1 = [11 11 11 11 11 11 11];
singleSpotRun2 = [9 8 7 6 5 3 1];

for kk=1:length(singleSpotRun1)
    allDifferences = compareSingleSpotsToDoubleSpot(singleSpotStatsFolder1, ...
        singleSpotRun1(kk), singleSpotStatsFolder2, singleSpotRun2(kk), ...
        doubleSpotStatsFolder, doubleSpotRun(kk));
    display(sprintf('Number of neurons missing spikes: %d; having too many: %d',...
        sum(allDifferences<0),sum(allDifferences>0)));
end

%% 2013-04-04-0

% % Dataset 0
% 
% singleSpotStatsFolder1 = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data001/statistics/activation';
% singleSpotStatsFolder2 = singleSpotStatsFolder1;
% doubleSpotStatsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data002/statistics/activation';
% 
% doubleSpotRun = [1:8];
% singleSpotRun1 = [1:8];
% singleSpotRun2 = [25:32];
% 
% for kk=1:length(singleSpotRun1)
%     allDifferences = compareSingleSpotsToDoubleSpot(singleSpotStatsFolder1, ...
%         singleSpotRun1(kk), singleSpotStatsFolder2, singleSpotRun2(kk), ...
%         doubleSpotStatsFolder, doubleSpotRun(kk));
% end

% % Dataset 1
% 
% singleSpotStatsFolder1 = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data003/statistics/activation';
% singleSpotStatsFolder2 = singleSpotStatsFolder1;
% doubleSpotStatsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/statistics/activation';
% 
% doubleSpotRun = [1:16];
% singleSpotRun1 = [1:8 1:8];
% singleSpotRun2 = [9:24];
% 
% for kk=1:length(singleSpotRun1)
%     allDifferences = compareSingleSpotsToDoubleSpot(singleSpotStatsFolder1, ...
%         singleSpotRun1(kk), singleSpotStatsFolder2, singleSpotRun2(kk), ...
%         doubleSpotStatsFolder, doubleSpotRun(kk));
% end

% % Dataset 2
% 
% singleSpotStatsFolder1 = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data003/statistics/activation';
% singleSpotStatsFolder2 = singleSpotStatsFolder1;
% doubleSpotStatsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data005/statistics/activation';
% 
% doubleSpotRun = [1:8];
% singleSpotRun1 = [9:16];
% singleSpotRun2 = [17:24];
% 
% for kk=1:length(singleSpotRun1)
%     allDifferences = compareSingleSpotsToDoubleSpot(singleSpotStatsFolder1, ...
%         singleSpotRun1(kk), singleSpotStatsFolder2, singleSpotRun2(kk), ...
%         doubleSpotStatsFolder, doubleSpotRun(kk));
% end