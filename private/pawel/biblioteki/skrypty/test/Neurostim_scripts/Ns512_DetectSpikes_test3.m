NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
%DataPath='E:\pawel\analysis\2009-11-27-0\data001\April2010\files'; 
DataPath= 'D:\Home\Rydygier\Neuro\files'; %'E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files'; %sciezka do danych RAW
AmplitudesVsChannelsAll_50us=[];
AmplitudesVsChannelsAll_100us=[];
RMSVsChannelsAll = [];

%usun szumiace kanaly:
BadChannels=[27:58 312];
Channels=[1:512];
Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=[1:64]; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]
Movies=[31:2:151]; %ktore amplitudy, domyslnie [1:2:151]

%Patterns=2;
%Movies=89;
%Channels=[442];
WritePathFigsGood='D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\50us_met12_good';
WritePathFigsBad='D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\50us_met12_bad';

for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,[1:512],NS_GlobalConstants);
    %[AmplitudesVsChannels,RMSVsChannels]=NS512_FindThresholdForLocalStimulation2PH(DataPath,WritePathFigs,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[6:75]); %7 to 76
    [AmplitudesVsChannels,RMSVsChannels]=NS512_FindThresholdForLocalStimulation3PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[6:75]); %7 to 76
    AmplitudesVsChannelsAll_50us=[AmplitudesVsChannelsAll_50us' AmplitudesVsChannels']';
    RMSVsChannelsAll = [RMSVsChannelsAll' RMSVsChannels']';
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