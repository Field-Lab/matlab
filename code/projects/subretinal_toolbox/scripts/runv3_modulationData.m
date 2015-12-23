clear;

neuronList = [];
psthStatsFolder = '/media/MEA_PROCESSED_8/2014-01-17-0/data/data001/statistics/all';
modulationStatsOutputFolder = '/media/MEA_PROCESSED_8/2014-01-17-0/data/data001/statistics/alternation_modulation';
alternationLogfilePath = '/media/MEA_PROCESSED_8/2014-01-17-0/logfiles/gratings_logfile001.log';
outputFolderFigures = '/media/MEA_PROCESSED_8/2014-01-17-0/data/data001/figures/alternation_modulation';

%% Computing statistics

stimNeuronList = computeAmplitudeModulation(psthStatsFolder, alternationLogfilePath, ...
    modulationStatsOutputFolder);

%% Plotting

allThresholds = plotAmplitudeModulation(modulationStatsOutputFolder, outputFolderFigures, ...
    'neuronList', stimNeuronList);

%% Getting the mean curve, normalized at 210um

% contentsModulationStatsFolder = dir(modulationStatsOutputFolder);
% neuronNames = struct('name','');
% nNeurons = 0;
% 
% for kk=1:length(contentsModulationStatsFolder)
%     if strfind(contentsModulationStatsFolder(kk).name,'.mat')
%         nNeurons = nNeurons + 1;
%         neuronNames(nNeurons).name = contentsModulationStatsFolder(kk).name;
%     end
% end
% 
% fh = figure; clf; set(fh,'color','white');
% hold on
% for kk=1:nNeurons
%     % Load the modulation data
%     load(fullfile(modulationStatsOutputFolder, neuronNames(kk).name));
%         
%     [gratings, ii] = sort(alternationData.gratings);
%     response = alternationData.response(ii);
%     response = response/response(gratings == 210);
%     
%     % Data
%     plot(gratings, response, '-^k', ...
%         'MarkerFaceColor',[0 0 0]);
% end

