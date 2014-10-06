FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='D:\Home\Rydygier\NEURO\files'; %sciezka do danych po preprocessingu
NS_GlobalConstants=NS_GenerateGlobalConstants(512);
%usun szumiace kanaly:
%BadChannels=[4 9 25 57];
Channels=[1:512];
%Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=[1:64]; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]

WritePathFigsGood='D:\Home\Pawel\analysis\retina\2009-11-27-0\AndrzejKoziec\2012-03-29-analysis\figures_good';
WritePathFigsBad='D:\Home\Pawel\analysis\retina\2009-11-27-0\AndrzejKoziec\2012-03-29-analysis\figures_bad';

% !!!!!!!!

MinDelay=6; %valid for 50 us pulses, for 100 us should be 6 

fid=fopen('D:\Home\Data\slices\2010-09-14-0\data002\data002000.bin','r');
Header=fread(fid,136);
fclose(fid);

cd D:\Home\Pawel\analysis\retina\2009-11-27-0\AndrzejKoziec\2012-03-29-analysis\data000;
nazwa = 'data000000.bin';
fid=fopen(nazwa,'a');
fwrite(fid, Header,'uint8');
fclose(fid);
      
% !!!!!!!!!!!
Movies=[1:150];
Movies=17;
ArrayID=500;
SamplesNumber=0;
SamplesOutput=0;
for i=1:1%length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,Channels,NS_GlobalConstants);
    [AmplitudesVsChannels,RMSVsChannels,SamplesOutput]=NS512_FindThresholdForLocalStimulation4AK_PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,[],[MinDelay:45],SamplesNumber); %7 to 76
    SamplesNumber=SamplesOutput;
end