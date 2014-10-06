FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='C:\Documents and Settings\Jojne\Pulpit\praktyki\data003'; %sciezka do danych po preprocessingu
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
%usun szumiace kanaly:
BadChannels=[4 9 25 57];
Channels=[1:64];
Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=Channels; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]

WritePathFigsGood='C:\Documents and Settings\Jojne\Pulpit\analiza\good';
WritePathFigsBad='C:\Documents and Settings\Jojne\Pulpit\analiza\bad';

% !!!!!!!!

MinDelay=6; %valid for 50 us pulses, for 100 us should be 6 

% !!!!!!!!!!!
Movies=[1:4];
ArrayID=1;
for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,[1:64],NS_GlobalConstants);
    [AmplitudesVsChannels,RMSVsChannels,Spikes_Array]=NS512_FindThresholdForLocalStimulation4PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,StimChannels,[MinDelay:40]); %7 to 76
 
end

