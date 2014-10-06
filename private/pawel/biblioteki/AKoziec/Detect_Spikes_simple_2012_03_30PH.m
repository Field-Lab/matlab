FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='D:\Home\Pawel\analysis\retina\2010-08-31-0\data009_preproc'; %sciezka do danych po preprocessingu
NS_GlobalConstants=NS_GenerateGlobalConstants(512);
Channels=[1:512];

[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,PatternsArray]=NS512_MoviePatterns('D:\Home\Data\retina\2010-08-31-0\movie009',NS_GlobalConstants);
Patterns=[65:320]; %ktore patterny (zestawy elektrod) chcemy przeanalizowac, domyslnie [1:64]

WritePathFigsGood='D:\Home\Pawel\analysis\retina\2010-08-31-0\2012-04-05-analysis\figures_good';
WritePathFigsBad='D:\Home\Pawel\analysis\retina\2010-08-31-0\2012-04-05-analysis\figures_bad';

% !!!!!!!!

MinDelay=6; %valid for 50 us pulses, for 100 us should be 6 

fid=fopen('D:\Home\Data\slices\2010-09-14-0\data002\data002000.bin','r');
Header=fread(fid,136);
fclose(fid);

cd D:\Home\Pawel\analysis\retina\2010-08-31-0\2012-04-05-analysis\data000;
nazwa = 'data000000.bin';
fid=fopen(nazwa,'a');
fwrite(fid, Header,'uint8');
fclose(fid);
      
% !!!!!!!!!!!
Movies=[1:270];
ArrayID=500;
SamplesNumber=0;
SamplesOutput=0;

for i=1:length(Patterns)
    PatternNumber=Patterns(i)
    Movies=find(PatternsArray(:,PatternNumber)==1)';
    %[StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,1,Channels,NS_GlobalConstants);
    [AmplitudesVsChannels,RMSVsChannels,SamplesOutput]=NS512_FindThresholdForLocalStimulation4AK_PH(DataPath,WritePathFigsGood,WritePathFigsBad,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,[],[MinDelay:45],SamplesNumber); %7 to 76
    SamplesNumber=SamplesOutput;
end