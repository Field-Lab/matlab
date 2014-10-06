FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='C:\Users\Jojne\Desktop\praktyki\praktyki\2009-11-27-0\files'; %sciezka do danych po preprocessingu
NS_GlobalConstants=NS_GenerateGlobalConstants(512);
%usun szumiace kanaly:
%BadChannels=[4 9 25 57];
Channels=[1:512];
%Channels=NS_RemoveBadChannels(Channels,BadChannels);

Patterns=[1:64]; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]

WritePathFigsGood='C:\Users\Jojne\Desktop\figures\good';
WritePathFigsBad='C:\Users\Jojne\Desktop\figures\bad';

% !!!!!!!!

MinDelay=6; %valid for 50 us pulses, for 100 us should be 6 

fid=fopen('C:\Users\Jojne\Desktop\praktyki\praktyki\2009-11-27-0\data\data001000.bin','r');
Header=fread(fid,136);
fclose(fid);
nazwa = 'data_050000.bin';
fid=fopen(nazwa,'a');
fwrite(fid, Header,'uint8');
fclose(fid);
      
% !!!!!!!!!!!
Movies=[1:151];
ArrayID=500;
SamplesNumber=0;
SamplesOutput=0;
for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,Channels,NS_GlobalConstants);
    [AmplitudesVsChannels,RMSVsChannels,SamplesOutput]=NS512_FindThresholdForLocalStimulation4PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,[],[MinDelay:40],SamplesNumber); %7 to 76
    SamplesNumber=SamplesOutput;
end

