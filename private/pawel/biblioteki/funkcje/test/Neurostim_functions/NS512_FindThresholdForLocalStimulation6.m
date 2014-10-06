function MovieNumber=NS512_FindThresholdForLocalStimulation5(DataPath,WritePathFigs,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
SpikesNumberThreshold=50;

ChannelsToRead=Channels;

for MovieNumber=Movies        
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);            
    DataTraces=DataTraces0(1:100,ChannelsToRead,Samples);
    
 %Wykrywanie spikow by Przemek
    %[TracesPrzemek, Artifact] = NS512_SubtractLocalArtifact(DataTraces, 10);
    %Events = FindSpikes_PR(TracesPrzemek, -25, 3, 3, 20);
 %Koniec wykrywania spikow by Przemek
    
 %Wykrywanie spikow metoda dwoch progow:
    [Events,Artifact]=NS512FindResponsesToLocalStimulation(DataTraces,25,3,10); %which traces include spikes that exceed threshold (independently for each channel)
    
    EventsPerChannel=sum(Events); %how many spikes on each electrode
    ChannelsWithSpikes=find(EventsPerChannel>SpikesNumberThreshold);    %Channels with minimum nimber of detected spikes
    
    for i=1:length(ChannelsWithSpikes)
        Channel=ChannelsToRead(ChannelsWithSpikes(i)); %na ktorym kanale wykryto spiki                               
        WaveformTypes=Events(:,ChannelsWithSpikes(i)); %0 - artifact only 1 - including spike
        artifacts= WaveformTypes==0; %find all the "artifact only" traces
        Artifact=mean(DataTraces(artifacts,:,:),1); % find EI of the artifact
        Traces=NS512_SubtractArtifact(DataTraces,Artifact); %Przebiegi po odjeciu artefaktow
        
        
        %Kawalek kodu by Przemek
        display('Elektrody sasiadujace z elektroda Channel');
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1)          %numery elektrody ze spikiem i elektrod sasiadujacych
        ChannelIndex = FindIndexes_PR(ChannelsToRead, ChannelsPlot);   %indeksy wskazujace numery elektrod w macierzy, 0 jesli dana elektroda jest wykluczona
       
        
        TracesSize = size(Traces);
        NumberOfAmplitudes = TracesSize(1,1);
        artifactsIndex = find(artifacts == 1);          %znajdywanie indeksow, w ktorych sa artefakty
        spikeIndex = find(artifacts == 0);              %znajdywanie indeksow, w ktorych sa spiki
       
        %Znajdywanie wektora wskaznikow dla spikow wybranych z maci. Traces
        ChannelTraces = Traces(:,ChannelIndex(1),:);    %Poszukiwanie ktora probka na elektrodzie na ktorej wykryto spike [ChannelIndex(1)]
                                                        %w danym movie jest pierwsza, ktora przekroczyla prog
        SCT=size(ChannelTraces);
        ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3));
        UniSpikesIndic = FindUnifiedSpikes_PR(ChannelTraces2D,-25) %Wektor numerow probek, ktore pierwsze przekraczaja prog w danym movie. 
        %Zwraca 0 jesli nie wystapilo przekroczenie progu
        
        TimeCorelSpikes = FindTimeCorellatedSpikes_PR(UniSpikesIndic, 2, 2); %Wektor zawierajacy 1 dla przebiegow zawartych w przedziale +/- k*sigma
        %0 dla przebiegow bez wykrytego spiku, 2 dla przebiegow poza
        %przedzialem +/-k*sigma od sredniego przesuniecia czasowego
        %wykrytych spikow na danej elektrodzie
        
        g1=find(TimeCorelSpikes==1);
        if length(g1)>0
            g2=min(UniSpikesIndic(g1));
            if g2>4
                UniSpikesIndic=UniSpikesIndic-4;
            else if g2>3 & g2<5
                %UniSpikesIndic=UniSpikesIndic-g2+1;
                UniSpikesIndic=UniSpikesIndic-3;
                else if g2>2 & g2<4
                        UniSpikesIndic=UniSpikesIndic-2;
                    else if g2>1 & g2<3
                            UniSpikesIndic=UniSpikesIndic-1;
                        end
                    end
                end
            end
        end

        N=7;
        figure(2)
        subplot(5,1,5);
        plot(sign(UniSpikesIndic),'bd-');
        %ChannelIndex = ChannelIndex(ChannelIndex~=0);
        ChannelTracesUUU = Traces(:,ChannelIndex,:);
        for i = 1:length(ChannelIndex)           
            if ChannelIndex(i) 
            ChannelTraces = Traces(:,ChannelIndex(i),:);              %przebiegi na elektrodzie o numerze Channel lub sasiadach
            ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3));
            subplot(5,N,i), h= plot(ChannelTraces2D');
            text(20,-80,num2str(ChannelsPlot(i)),'Fontsize',16);
            set(h(artifactsIndex),'Color','Black');
            set(h(spikeIndex),'Color','Red');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
            
            subplot(5,N,i+N), ddd = plot(ChannelTraces2D'); %WithoutArtefacts');
            set(ddd(find(TimeCorelSpikes==0)),'Color','Black');
            set(ddd(find(TimeCorelSpikes==1)),'Color','Red');
            set(ddd(find(TimeCorelSpikes==2)),'Color','Blue');
            axis([0 40 -100 50]);
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
            grid on;
            
            
            
            
%            UniSpikesIndic = UniSpikesIndicator(find(UniSpikesIndic & (TimeCorelSpikes == 1)));
            [SpikeUnif, MeanSpike] = SpikeUnif_PR(ChannelTraces2D, UniSpikesIndic, TimeCorelSpikes);
            subplot(5,N,i+2*N), f = plot(SpikeUnif');
            set(f,'Color','Red');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
            subplot(5,N,i+3*N), g = plot(MeanSpike');
            set(g,'Color','Blue');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
%             ChannelTraces = Traces(:,ChannelIndex(1),:);
%             ChannelTraces2D = reshape(ChannelTraces,NumberOfAmplitudes,40);
%             display('Rozmiar macierzy ChannelTraces2D');
%             [SpikeUnif, MeanSpike] = SpikeUnif_PR(ChannelTraces2D);
%             subplot(2,7,8), g = plot(SpikeUnif');
%             subplot(2,7,9), f = plot(MeanSpike');
%             set(g,'Color','Red');
%             axis([0 40 -100 50]);
%             grid on;
            end
        end
        
        
        hc=gcf;
        FullName=[WritePathFigs '\' 'local_' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber) '_el' num2str(Channel)];
        set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[16 9]);
        set(hc,'PaperPosition',[0 0 16 9]);  
        print(hc, '-dtiff', '-r120', FullName);
        %Koniec kawalka kodu by Przemek
        
        %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(Traces,ChannelsToRead,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Channel);                
        display('rozmiary macierzy ChannelTraces, ChannelsPlot, UniSpikesIndicator');
            size(ChannelTracesUUU)
            ChannelsPlot %= ChannelsPlot(ChannelsPlot~=0)';
            size(ChannelsPlot)
            size(UniSpikesIndic)
        
        %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks_PR(ChannelTracesUUU,ChannelsPlot,TimeCorelSpikes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Channel);
        h=gcf;
        FullName=[WritePathFigs '\' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber) '_el' num2str(Channel)];            
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        %print(h, '-dtiff', '-r120', FullName);
    end
    ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,ChannelsToRead(ChannelsWithSpikes));
end
%if a<51
%    MovieNumber=0;
%end