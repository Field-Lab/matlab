StimulatedElectrodes=StimulatedElectrodesMouse;
s1=size(StimulatedElectrodes);

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
for i=10:length(PatternsToAnalyze)
    clf
    StimulatedChannel=PatternsToAnalyze(i);
    ChannelData=StimulatedElectrodes(:,StimulatedChannel,:);
    qw(StimulatedChannel)=max(max(max(ChannelData)));
    
    ChannelDataNew=reshape(ChannelData,s1(1),s1(3));
    ChannelsToExclude=electrodeMap.getAdjacentsTo(StimulatedChannel,1)';
    ChannelDataNew(:,ChannelsToExclude)=0; %ignore spikes on the stimulated electrode
    
    qw2(StimulatedChannel)=max(max(max(ChannelData)));
    
    for Amplitude=1:25
        subplot(5,5,Amplitude);        
        h1=plot(X(StimulatedChannel),Y(StimulatedChannel),'rd');
        set(h1,'MarkerSize',6);
        set(h1,'MarkerFaceColor','r');
        hold on;
        Amplitude
        for RecEl=1:512
            h1=plot(X(RecEl),Y(RecEl),'bo');            
            set(h1,'MarkerSize',4);
            %ColorIndex=ChannelDataNew(Amplitude,RecEl)*5/50+1;
            set(h1,'MarkerEdgeColor',MapaKolorow(ChannelDataNew(Amplitude,RecEl)+1,:));
            set(h1,'MarkerFaceColor',MapaKolorow(ChannelDataNew(Amplitude,RecEl)+1,:));            
            %set(h1,'MarkerSize',round(ChannelDataNew(Amplitude,RecEl)*5/50)+1);                     
        end   
        axis([-1000 1000 -500 500])
        h=gca;        
        set(h,'Color',[0 0 0]);
        set(h,'XTickLabel','');
        set(h,'YTickLabel',''); 
        set(h,'Box','on');
    end    
    FullName=[FigurePath '\Pattern' num2str(StimulatedChannel) '.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    print(h, '-dtiff', '-r120', FullName);
end

