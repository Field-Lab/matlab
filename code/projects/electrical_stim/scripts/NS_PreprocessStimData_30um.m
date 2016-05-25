function NS_PreprocessStimData_30um(pathToRawData,WritePath)
% Function to preprocess electrical stimulation data acquired on the 30 um
% stimulation board, SB512-3. This preprocessing step uses a transform
% function to map the electrodes acquired using the 512-electrode LabVIEW 
% code to the 519 electrode map. 2016-05-23. Hopefully it will one day be
% updated to avoid this transform.
% Inputs: strings to the raw data and save paths, for example: 
% pathToRawData = '/Volumes/Lab/Transfer/30umtest/data000';
% WritePath='/Volumes/Lab/Transfer/30umtest/processed/data000-2/';

% Generate structure containing global constants
ChipAddresses=24:31;
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

% Set array ID for the 519/30um stim board
ArrayID=1502;

% Preprocess the data. 
WritePathFigs = WritePath;
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
ChannelsToRead = 1:512;
makePlots = 0; % Choose whether or not to save plots. 
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=...
    NS_PreprocessDataNew512v6main519(pathToRawData,WritePath,WritePathFigs,...
ArrayID,FigureProperties,NS_GlobalConstants,makePlots,[1:512],ChannelsToRead,0);