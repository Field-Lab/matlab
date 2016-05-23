NS_GlobalConstants=NS_GenerateGlobalConstants(500);
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatterns,Patterns]=NS512_MoviePatterns('J:\data\2013-12-12-3-PH\movie001',NS_GlobalConstants);

StimulatedElectrodes=StimulatedElectrodesMouse;
s1=size(StimulatedElectrodes);
break
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);        
end

FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2013-12-12-3-PH-2014-08-16';

MapaKolorow=jet;
clear qw;
PatternsToAnalyze=unique(AllPatterns)';

figure(1)
for i=17%:length(PatternsToAnalyze)
    clf
    StimulatedChannel=PatternsToAnalyze(i);
    ChannelData=StimulatedElectrodes(:,StimulatedChannel,:);
    qw(StimulatedChannel)=max(max(max(ChannelData)));
    
    ChannelDataNew=reshape(ChannelData,s1(1),s1(3));
    ChannelsToExclude=electrodeMap.getAdjacentsTo(StimulatedChannel,1)';
    ChannelDataNew(:,ChannelsToExclude)=0; %ignore spikes on the stimulated electrode
    
    qw2(StimulatedChannel)=max(max(max(ChannelData)));
    
    for Amplitude=2:25
        Amplitude;
        %subplot(4,6,Amplitude-1);
        SubplotRow=ceil((Amplitude-1)/6)
        SubplotColumn=Amplitude-1-(SubplotRow*6)+6-1;
        subplot('position',[SubplotColumn*0.165+0.008 (4-SubplotRow)*0.24+0.03 0.16 0.22]);
        %h1=plot(X(StimulatedChannel),Y(StimulatedChannel),'rd');
        %set(h1,'MarkerSize',6);
        %set(h1,'MarkerFaceColor','r');
        %hold on;
        
        ChannelDataForAmplitude=ChannelDataNew(Amplitude,:);
        ConvertedData=ChannelDataToColors(ChannelDataForAmplitude,X,Y);
        image(ConvertedData');
        colormap(jet(50))
        hold on
        h10=plot(X(StimulatedChannel)/30+33,17-round((Y(StimulatedChannel)+490)/60),'go');
        set(h10,'MarkerSize',8)
        %axis([-1000 1000 -500 500])
        h=gca;        
        %set(h,'Color',[0 0 0]);
        set(h,'XTickLabel','');
        set(h,'YTickLabel',''); 
        %et(h,'Box','on');
    end    
    FullName=[FigurePath '\Pattern' num2str(StimulatedChannel) 'v2.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 5.5]);
    set(h,'PaperPosition',[0 0 16 5.5]); 
    %print(h, '-dtiff', '-r120', FullName);
end

