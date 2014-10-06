NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ArrayID=500;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
Radius=1;

FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-100 50],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
DataPath='I:\analysis\retina\2009-11-27-0\data001';

%PatternNumber=20;
%MovieNumber=57;

Patterns=[1:64];

WritePathFigs='I:\analysis\retina\2009-11-27-0\data001\figures2\50us';
Movies=[35:2:151];
%Channels=[1:12];
BadChannels=[312];

for i=1:length(Patterns)
    Pattern=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,Pattern,1,[1:512],NS_GlobalConstants);
    for j=8 %1:length(StimChannels)
        CenterChannel=StimChannels(j)
        CenterChannels(Pattern,j)=CenterChannel;
        Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';        
        Channels=[385:512];
        Channels=NS_RemoveBadChannels(Channels,BadChannels);
        MarkedChannels=[CenterChannel-64 CenterChannel];
        MovieNumber=NS512_FindThresholdForLocalStimulation(DataPath,Pattern,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels);
        Amplitudes(Pattern,j)=MovieNumber;
        if MovieNumber~=0
            FullName=[WritePathFigs '\' 'el' num2str(CenterChannel) '_and' num2str(CenterChannel-64) '_p' num2str(Pattern) '_m' num2str(MovieNumber)];
            hj=gcf;
            set(hj, 'PaperOrientation', 'portrait');
            print(hj, '-dtiff', FullName);
        end
    end
end

break;
WritePathFigs='D:\analysis\retina\2009-11-27-0\data001\figures\100us';
Movies=[2:2:150];
%Channels=[1:12];
BadChannels=[312];

for i=1:length(Patterns)
    PatternNumber=Patterns(i);
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,Patterns(i),1,[1:512],NS_GlobalConstants);
    for j=1:length(StimChannels)
        CenterChannel=StimChannels(j)
        CenterChannels(i,j)=CenterChannel;
        Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';        
        Channels=NS_RemoveBadChannels(Channels,BadChannels);
        MovieNumber=NS512_FindThresholdForLocalStimulation(DataPath,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants);
        Amplitudes(i,j)=MovieNumber;
        if MovieNumber~=0
            FullName=[WritePathFigs '\' 'el' num2str(CenterChannel)  '_p' num2str(Patterns(i)) '_m' num2str(MovieNumber)];
            hj=gcf;
            set(hj, 'PaperOrientation', 'portrait');
            print(hj, '-dtiff', FullName);
        end
    end
end