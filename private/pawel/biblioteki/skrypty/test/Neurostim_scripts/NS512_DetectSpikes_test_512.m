NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArrayID=1;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
%DataPath='E:\pawel\analysis\2009-11-27-0\data001\April2010\files'; 
DataPath= 'C:\home\pawel\2010\analysis\retina_61\2010-08-20-0\nowe\data007'; %sciezka do danych po preprocessingu
AmplitudesVsChannelsAll_50us=[];
AmplitudesVsChannelsAll_100us=[];
RMSVsChannelsAll = [];

%usun szumiace kanaly:
BadChannels=[4 9 25 57];
Channels=[1:64];
Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=Channels; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]
Movies=[1:26]; %ktore amplitudy, domyslnie [1:2:151]

%Patterns=2;
%Movies=89;
%Channels=[442];
WritePathFigsGood='C:\home\pawel\2010\analysis\retina_61\2010-08-20-0\nowe\50us_good';
WritePathFigsBad='C:\home\pawel\2010\analysis\retina_61\2010-08-20-0\nowe\50us_bad';

%Channels=[2 10 17 20 27];
%Channels=[1 2 7 10 15 16 17 19 23 27 28 63];
%Channels=[1 28 63];
Thresholds=[30 30 60 30 40];

for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,[1:64],NS_GlobalConstants);
    %[AmplitudesVsChannels,RMSVsChannels]=NS512_FindThresholdForLocalStimulation2PH(DataPath,WritePathFigs,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[6:75]); %7 to 76
    [AmplitudesVsChannels,RMSVsChannels]=NS512_FindThresholdForLocalStimulation4PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[6:50]); %7 to 76
    %AmplitudesVsChannelsAll_50us=[AmplitudesVsChannelsAll_50us' AmplitudesVsChannels']';
    %RMSVsChannelsAll = [RMSVsChannelsAll' RMSVsChannels']';
end
break;
cd D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures;
save('AmplitudesVsChannels_50us.mat','AmplitudesVsChannelsAll_50us');
%f=fopen('AmplitudesVsChannels_50us','w+');
%fwrite(f,AmplitudesVsChannelsAll);
%fclose(f);
f=fopen('RMSVSChannels_50us','w+');
fwrite(f,RMSVsChannelsAll);
fclose(f);