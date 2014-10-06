function MovieNumber=NS512_FindThresholdForLocalStimulation7a(DataPath,WritePathFigs,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
SpikesNumberThreshold=50;

ChannelsToRead=Channels;

for MovieNumber=Movies      
  %Wczytanie danych
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);            
    DataTracesFull=DataTraces0(1:100,[1:512],Samples);
    
    DataTraces=DataTraces0(1:100,ChannelsToRead,Samples); % wyluskaj dane z "dobrych" elektrod
 
 % * * *
 %[TracesPrzemek, Artifact] = NS512_SubtractLocalArtifact(DataTraces, 10);    %TracesPrzemek - surowe dane bez artefaktow
    [TracesWithoutArtifact, Artifact] = NS512_SubtractLocalArtifact(DataTracesFull, 10);
    TracesWithoutArtifactLimited=TracesWithoutArtifact(:,ChannelsToRead,:);  
    
    STWA=size(TracesWithoutArtifact)
    Events=zeros(STWA(1),STWA(2));
    Events0 = FindSpikes_PR(TracesWithoutArtifactLimited, -25, 3, 3, 20); %Events - macierz dwuwymiarowa (liczba elektrod razy liczba przebiegow)
    size(Events0)
    Events(:,ChannelsToRead)=Events0; % Events - rozmiar 512 x 100; Events0 - rozmiar 479 x 100
            
    EventsPerChannel=sum(Events); %how many spikes on each electrode
    ChannelsWithSpikes=find(EventsPerChannel>SpikesNumberThreshold) % ktore elektrody zarejestrowaly wiecej niz 50 spikow
 
    for i=1:length(ChannelsWithSpikes) % dla wszystkich elektrod z wiecej niz 50 spikami...
        Channel=ChannelsWithSpikes(i);
        WaveformTypes=Events(:,ChannelsWithSpikes(i)); %0 - artifact only 1 - including spike
        artifacts= WaveformTypes==0; %find all the "artifact only" traces
        Artifact=mean(DataTraces(artifacts,:,:),1); % find EI of the artifact
        Traces=NS512_SubtractArtifact(DataTraces,Artifact); %Przebiegi po odjeciu artefaktow
                
        %Kawalek kodu by Przemek
        display('Elektrody sasiadujace z elektroda Channel');
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1)
        ChannelIndex = FindIndexes_PR(Channels, ChannelsPlot)    %indeksy wskazujace numery elektrod
        clf;
        
        TracesSize = size(Traces);
        NumberOfAmplitudes = TracesSize(1,1);
        artifactsIndex = find(artifacts == 1);
        spikeIndex = find(artifacts == 0);
       
        %Znajdywanie wektora wskaznikow dla spikow wybranych z maci. Traces
        ChannelTraces = Traces(:,ChannelIndex(1),:);
        SCT=size(ChannelTraces);
        ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3));
        UniSpikesIndic = FindUnifiedSpikes_PR(ChannelTraces2D,-25); %Wektor wskaznikow (numery probek dla ktorych nastapilo przekroczenie progu)
        TimeCorelSpikes = FindTimeCorellatedSpikes_PR(UniSpikesIndic, 2, 2); %Wektor zawierajacy 1 dla przebiegow zawartych w przedziale +/- k*sigma               
       
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
            
            %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(ChannelTraces,ChannelsWithSpikes,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Channel);
            
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