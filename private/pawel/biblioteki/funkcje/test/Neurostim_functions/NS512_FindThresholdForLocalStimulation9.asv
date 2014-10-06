function [AmplitudesVsChannels, RMSVsChannelsAll ] = NS512_FindThresholdForLocalStimulation9(DataPath,WritePathFigs,PatternNumber,Movies,Channels,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Samples);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
SpikesNumberThreshold=50;

ChannelsToRead=Channels;
AmplitudesVsChannels=[];
RMSVsChannelsAll = [];

for MovieNumber=Movies      
  %Wczytanie danych
  %display(MovieNumber)
  [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);            
    DataTracesFull=DataTraces0(1:100,[1:512],Samples);
    
    Test1 = [];
    Test2 = [];
    [Test1, Test2] = FindTracesClasses_PR(DataTracesFull, 10, 50);
    
    gkgjhk=size(Test1)
    gkgk=size(Test2)
    figure(103);
   
    
    DataTraces=DataTraces0(1:100,ChannelsToRead,Samples); % wyluskaj dane z "dobrych" elektrod
 
    [TracesWithoutArtifact, Artifact] = NS512_SubtractLocalArtifact(DataTracesFull, 10);
    TracesWithoutArtifactLimited=TracesWithoutArtifact(:,ChannelsToRead,:);  
    
    STWA=size(TracesWithoutArtifact);
    Events=zeros(STWA(1),STWA(2));
    [Events0, SpikesBeyondTimeWindow] = FindSpikes_PR(TracesWithoutArtifactLimited, -25, 3, 3, 20); %Events - macierz dwuwymiarowa (liczba elektrod razy liczba przebiegow)
    size(Events0)

    Events(:,ChannelsToRead)=Events0; % Events - rozmiar 512 x 100; Events0 - rozmiar 479 x 100
            
    EventsPerChannel=sum(Events); %how many spikes on each electrode, EventsPerChannel - rozmiar 1 x 512
    
    EventsBeyondTimeWindow = zeros(STWA(1),STWA(2));
    EventsBeyondTimeWindow(:,ChannelsToRead) = SpikesBeyondTimeWindow;
%     EventsBeyondTimeWindow = sum(EventsBeyondTimeWindow==2);    %znajduje, ile wystapilo przekroczen progu zbyt wczesnie lub zbyt pozno
     %na danej elektrodzie, EventsBeyondTimeWindow - macierz o wymiarach 1 x 512 
    
    ChannelsWithSpikes=find(EventsPerChannel > SpikesNumberThreshold); % ktore elektrody zarejestrowaly wiecej niz 50 spikow
        

    for i=1:length(ChannelsWithSpikes) % dla wszystkich elektrod z wiecej niz 50 spikami...
        Channel=ChannelsWithSpikes(i);
        WaveformTypes=Events(:,ChannelsWithSpikes(i)); %0 - artifact only 1 - including spike
        artifacts= WaveformTypes==0; %find all the "artifact only" traces
        
        ArtifactLimited=mean(DataTraces(artifacts,:,:),1); % find EI of the artifact (479 kanalow)
        Traces=NS512_SubtractArtifact(DataTraces,ArtifactLimited); %Przebiegi po odjeciu artefaktow (479 kanalow)
        
        Artifact=mean(DataTracesFull(artifacts,:,:),1); % find EI of the artifact (512 kanalow)
        TracesFull=NS512_SubtractArtifact(DataTracesFull,Artifact); %Przebiegi po odjeciu artefaktow (512 kanalow)
                
        %Kawalek kodu by Przemek
        ChannelsPlot = electrodeMap.getAdjacentsTo(Channel,1); % Elektrody sasiadujace z elektroda Channel (ta z wiecej niz 50 spikami)
        %ChannelIndex = FindIndexes_PR(Channels, ChannelsPlot);    %indeksy wskazujace numery elektrod
        ChannelIndex=ChannelsPlot;      %Numery elektrod sasiadujacych z elektroda Channel
        clf;
        
        TracesSize = size(Traces);
        NumberOfAmplitudes = TracesSize(1,1);        
        
        artifactsIndex = find(WaveformTypes == 0);      %Wskazniki, ktore przebiegi sa artefaktami
        spikeIndex = find(WaveformTypes == 1);          %Wskazniki, ktore przebiegi sa spikami
       
        %Znajdywanie wektora wskaznikow dla spikow wybranych z maci. Traces
        ChannelTraces = TracesFull(:,ChannelIndex(1),:);    %Wszystkie 100 przebiegow dla wybranej elektrody
        SCT=size(ChannelTraces);
        ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3)); %100 przebiegow przeksztalcone w macierz 2D
        UniSpikesIndic = FindUnifiedSpikes_PR(ChannelTraces2D,-25); %Wektor wskaznikow (numery probek dla ktorych nastapilo przekroczenie progu)       
        TimeCorelSpikes = FindTimeCorellatedSpikes_PR(UniSpikesIndic, 1.5, 1.5); %Wektor zawierajacy 1 dla przebiegow zawartych w przedziale +/- k*sigma
        
        %Obliczanie liczby spikow na danej elektrodzie wykrytych poza oknem czasowym
        %zdefiniowanym w funkcji FindSpikes_PR
        TracesNumberBeyondTimeWindow = (EventsBeyondTimeWindow(:,ChannelIndex(1)).*(~TimeCorelSpikes))+EventsBeyondTimeWindow(:,ChannelIndex(1));
        TracesNumberBeyondTimeWindow = TracesNumberBeyondTimeWindow>0;
        TracesNumberBeyondTimeWindow = sum(TracesNumberBeyondTimeWindow);
        
        if EventsPerChannel(Channel)- TracesNumberBeyondTimeWindow < 40
            display(strcat('elektroda ', num2str(ChannelIndex(1)), ' wylaczona'));
            continue
        end
        
        g1=find(TimeCorelSpikes==1);         
        if length(g1)>0 % jesli sa jakies spiki z niewielkim jitterem
            g2=min(UniSpikesIndic(g1));
            if g2>4
                UniSpikesIndic(g1)=UniSpikesIndic(g1)-4;
            else
                UniSpikesIndic(g1)=UniSpikesIndic(g1)-g2+1;
            end
        end     
                
        N=7;
%{
        figure(2)
        subplot(5,2,9);
        plot(sign(UniSpikesIndic),'bd-');
 %}
        for i = 1:length(ChannelIndex)           
            if ChannelIndex(i) 
            ChannelTraces = TracesFull(:,ChannelIndex(i),:); %przebiegi na elektrodzie o numerze Channel lub sasiadach
            ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3));
%{
            subplot(5,N,i), h= plot(ChannelTraces2D');
            text(20,-80,num2str(ChannelsPlot(i)),'Fontsize',16);
            set(h(artifactsIndex),'Color','Black');
            set(h(spikeIndex),'Color','Red');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
           
            subplot(5,N,i+N), ddd = plot(ChannelTraces2D'); %WithoutArtifacts;
            set(ddd(find(TimeCorelSpikes==0)),'Color','Black');
            set(ddd(find(TimeCorelSpikes==1)),'Color','Red');
            set(ddd(find(TimeCorelSpikes==2)),'Color','Blue');
            axis([0 40 -100 50]);
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
            grid on;
%}            
            %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(ChannelTraces,ChannelsWithSpikes,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Channel);
            
            [SpikeUnif, MeanSpike] = SpikeUnif_PR(ChannelTraces2D, UniSpikesIndic, TimeCorelSpikes); %Wyznaczanie macierzy spikow 
                                                                                                %uporzadkowanych czasowo oraz spiku usrednionego

        %Obliczanie sredniej z RMS z kazdego przebiegu
            for m=1:size(SpikeUnif, 1)      %Obliczanie roznicy przebiegow od spiku usrednionego
            TracesMinusMeanSpike(m,:) = SpikeUnif(m,:) - MeanSpike;
            TracesRMS(m) = std(TracesMinusMeanSpike(m,:));  %RMS wszystkich przebiegow zakwalifikowanyhc jako spikei
            end
            
            MeanTracesRMS = mean(TracesRMS);    %Srednia wartosc odchylen od przebiegu usrednionego
            MeanSpikePeakToPeak = max(MeanSpike) - min(MeanSpike);  %Amplituda PTP sredniego spiku
            RMSPercentOfMeanSpike = 100*MeanTracesRMS/MeanSpikePeakToPeak;  %roznica procentowa do umieszczenia na rysunku
            
            if i==1                     %Dla pierwszej rysowanej elektrody utworz zmienna RMSVsChannels do zapisania do pliku
            RMSVsChannels = [PatternNumber MovieNumber ChannelIndex(1) MeanSpikePeakToPeak MeanTracesRMS];
            end
 %{           
            subplot(5,N,i+2*N), f = plot(SpikeUnif');
            text(20,-80,sprintf('%0.3g',MeanTracesRMS),'Fontsize',16);
            set(f,'Color','Red');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
            
            subplot(5,N,i+3*N), g = plot(MeanSpike');
            text(20,-80,strcat(sprintf('%0.4g',RMSPercentOfMeanSpike),'%'),'Fontsize',16);
            set(g,'Color','Blue');
            axis([0 40 -100 50]);
            grid on;
            h23=gca;
            set(h23, 'XTick', [0:5:40]);
            set(h23, 'YTick', [-100:20:40]);
%}            
            end
        end
        
     %Tworzenie macierzy w celu zapisania do pliku
        EI=NS512_EI_Corrected_Timings(TracesFull,[1:512],1,UniSpikesIndic,TimeCorelSpikes,30);
        AmplitudeVsChannel=max(EI')-min(EI');
        AmplitudesVsChannels=[AmplitudesVsChannels' AmplitudeVsChannel']';
        
        RMSVsChannelsAll = [RMSVsChannelsAll; RMSVsChannels];
     
     %rysowanie zarejestrowanych amplitud w funkcji kanalow
%{
        subplot(5,2,10);
        plot(AmplitudeVsChannel,'bd-');
        grid on;
        h=gca;
        set(h,'XLim',[-10 520]);
        
        hc=gcf;
        FullName=[WritePathFigs '\' 'local_' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber) '_el' num2str(Channel)];
        set(hc,'PaperUnits','inches');
        set(hc,'PaperSize',[16 9]);
        set(hc,'PaperPosition',[0 0 16 9]);  
        print(hc, '-dtiff', '-r120', FullName);
        %Koniec kawalka kodu by Przemek
%}
        %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(Traces,ChannelsToRead,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants,MarkedChannels,Channel);                        
        h=gcf;
        FullName=[WritePathFigs '\' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber) '_el' num2str(Channel)];            
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        %print(h, '-dtiff', '-r120', FullName);
    end
    %ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,ChannelsToRead(ChannelsWithSpikes));
    ChannelsToRead=NS_RemoveBadChannels(ChannelsToRead,ChannelsWithSpikes);
end
%if a<51
%    MovieNumber=0;
%end