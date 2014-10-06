clc
clear all
close all

set(0, 'defaultfigurecolor', 'w')
%%% ELEKTRODY %%%%%% 
ArrayChannelRead=[1:3,5:8,10:24,26:56,58:64];   %%recznie wygenerowana tablica odczytu

%%%%%% ERROR 
CutOff =4;    %poziom odciecia

%%%% PATHS %%%%%
javaaddpath 'C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\Vision.jar'      %%% !!!!!! do zmiany przy innym komp
DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023';   %%% !!!!!! do zmiany przy innym komp
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);     %define the electrode map - must be different than 1 for silicon probes
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;

LicznikParElektrod=1;
Radius=1;  % odleg³oœæ 60 mikronów (2-120,3-180)
OverTreshold(1:25)=0;
OdKtorego=OverTreshold(1:25);
ColorScale=[0,0,0.5977;0,0,0.6953;0,0,0.7969;0,0,0.8945;0,0,0.9961;0,0.1106,0.9961;0,0.1992,0.9961;0,0.3008,0.9961;0,0.3984,0.9961;0,0.5,0.9961;0,0.5977,0.9961;0,0.6953,0.9961;0,0.7969,0.9961;0,0.8945,0.9961;0,0.9961,0.9961;0.1016,0.9961,0.8945;0.1992,0.9961,0.7969;0.3008,0.9961,0.6953;0.3984,0.9961,0.5977;0.5,0.9961,0.5;0.5977,0.9961,0.3984;0.3008,0.9961,0.6953;0.7969,0.9961,0.1992;0.8945,0.9961,0.1016;0.9961,0.9961,0;0.9961,0.8945,0;0.9961,0.7969,0;       0.9961,0.6953,0;0.9961,0.7969,0;0.9961,0.5977,0;0.9961,0.5,0;0.9961,0.3984,0;0.9961,0.3008,0;0.9961,0.1992,0;0.9961,0.1016,0;0.9961,0,0;0.8945,0,0;0.7969,0,0;0.5977,0,0;0.5,0,0];

%%%%% ILOSC POMIAROW
NoOfHistograms=25;
NoOfMovies=25;
ElectrodeOnFig=63;
ProblematycznePary=zeros(2,50);
ProblematyczneParyCounter=1;
WhenTracesUnderCutOff=8;   %% od ktorej probki ustawiamy gdy przeieg ponizej cutoff

%%%%% ANALIZA MOVIE
MovieDoAnalizy=11;   % na wykresie coordinates wyœwietlone polaczenia w przypadku osi¹gniêcia wybranej wartoœci
OdKtorego_Analisys=10;  % od ktorego przypadku uznajemy wartosc OdKtorego za zbyt duza

%%%% COORDINATES  - rysowanie matrycy elektrod

PlotCoordinates(ElectrodeOnFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
      
      
for ElectrodeArrayCounter=1:1  % max(size(ArrayChannelRead))   %%60 elektrod (bez 4,9,25,57)    
   
    CenterChannel = ArrayChannelRead(ElectrodeArrayCounter);
    ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';  
    ChannelsPlot_BezWlasnej=ChannelsPlot(2:max(size(ChannelsPlot))); %%wersja bez dzialania elektrody na elektrode    
           
    for j=1:max(size(ChannelsPlot_BezWlasnej))     %%%%% pierwotnie: :max(size(ChannelsPlot))

        ChannelRead=ChannelsPlot_BezWlasnej(j);         
        MovieNumber = 1;                

        [DataTraces]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
        SDT=size(DataTraces);
        DataTraces=DataTraces(1:23,:,1:40);
        ClusterWithArtifacts=[1:SDT(1)]';             
                
        ArtifactsMatrix=size(ClusterWithArtifacts);
        Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
        Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
     
        FirstArtifact=Przebieg1;
            
        NS_GlobalConstans=NS_GenerateGlobalConstants(61);
        [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
        FirstMovieAmplitude=Amplitudes;
        %%%% SECOND ARTIFACT

        MovieNumber = 2;
        
        [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
        SDT=size(DataTraces);
        DataTraces=DataTraces(1:23,:,1:40);
        ClusterWithArtifacts=[1:SDT(1)]';             
                                
        ArtifactsMatrix=size(ClusterWithArtifacts);
        Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
        Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)'); 
             
        SecondArtifact=Przebieg2;       
          
        NS_GlobalConstans=NS_GenerateGlobalConstants(61);
        [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
        SecondMovieAmplitude=Amplitudes;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        C2=(FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude);
        C1=FirstArtifact-((FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude)).*FirstMovieAmplitude;
        
        for MovieNumberCounter=1:NoOfMovies    %%25 movies              
            
            MovieNumber = MovieNumberCounter;            
                
            [DataTraces]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
            SDT=size(DataTraces);
            DataTraces=DataTraces(1:23,:,1:40);
            ClusterWithArtifacts=[1:SDT(1)]';

            ArtifactsMatrix=size(ClusterWithArtifacts);
            Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
            Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');

            Avg_Przebieg1=Przebieg1;

            MovieNumber = MovieNumberCounter+1; 

            [DataTraces]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
            SDT=size(DataTraces);
            DataTraces=DataTraces(1:23,:,1:40);
            ClusterWithArtifacts=[1:SDT(1)]';

            ArtifactsMatrix=size(ClusterWithArtifacts);
            Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
            Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');

            Avg_Przebieg2=Przebieg2;

             C2=(FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude);
             C1=FirstArtifact-((FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude)).*FirstMovieAmplitude;
            
           if( MovieNumberCounter < 4)              
                                     
                MovieNumber =MovieNumberCounter;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
               [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

                MovieNumber =MovieNumberCounter+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;            

                               
            else                

                MovieNumber=MovieNumberCounter-2;

                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                
                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';                                
                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
                Avg_Przebieg1_fixed=mean(Przebieg1');           

                FirstArtifact=Przebieg1;
                
                MovieNumber=MovieNumberCounter-1;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;   

                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';              
                                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');
                Avg_Przebieg2_fixed=mean(Przebieg2');              

                SecondArtifact=Przebieg2;
                
                C2=(FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude);
                C1=FirstArtifact-((FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude)).*FirstMovieAmplitude;

                MovieNumber =MovieNumberCounter;
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;  

                MovieNumber =MovieNumberCounter+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;           
                               
            end      
                Wykres=(Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt);
                Wykres=Wykres';
                Wykres(:,1:7)=0;   %% zerowanie pierwszych 7 probek
        
                WykresAvg=mean(Wykres);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
               
                if (find((abs(WykresAvg))>CutOff))
                    OverTresholdMatrix=find((abs(WykresAvg))> CutOff );                 
                    OdKtorego(MovieNumberCounter)=max(OverTresholdMatrix); %% od ktorej probki wykres jest ponizej cutoff
                else
                     
                    OdKtorego(MovieNumberCounter)=WhenTracesUnderCutOff;   %% jezeli find nic nie wykryje automatycznie ustawiane jest od 8 (pierwszej)probki
                end
                
                Elektroda_Stymulujaca=ArrayChannelRead(ElectrodeArrayCounter);   %% polaczenia elektrod na wykresie coordinates
                Elektroda_Odczytowa=ChannelRead;

                 if (MovieNumber == MovieDoAnalizy && OdKtorego(MovieNumberCounter) > OdKtorego_Analisys);   
                      
                    ProblematycznePary(:,ProblematyczneParyCounter)=[Elektroda_Stymulujaca,Elektroda_Odczytowa];
                    
                    figure (1)
                    hold on
                    X_Coords=[electrodeMap.getXPosition(Elektroda_Stymulujaca),electrodeMap.getXPosition(Elektroda_Odczytowa)];
                    Y_Coords=[electrodeMap.getYPosition(Elektroda_Stymulujaca),electrodeMap.getYPosition(Elektroda_Odczytowa)];
                    hold on                   
                    
                    plot(X_Coords,Y_Coords,'Color',ColorScale(OdKtorego(MovieNumberCounter),:))
                    text(electrodeMap.getXPosition(Elektroda_Odczytowa),electrodeMap.getYPosition(Elektroda_Odczytowa), '+','FontSize',15,'VerticalAlignment','middle','HorizontalAlignment','center', 'Color',ColorScale(OdKtorego(MovieNumberCounter),:))  %% '+' pojawia sie na el odczytowej    
                    colorbar('YTickLabel',{'5','10','15','20','25','30','35','40'})
                    hold on  
                    

                    figure (4)
                    subplot(5,2,ProblematyczneParyCounter)    %% mozna zmienic wartosc 5
                    title(['Stim. el: ',num2str(ProblematycznePary(1,ProblematyczneParyCounter)),'\newlineRead el. : ', num2str(ProblematycznePary(2,ProblematyczneParyCounter))])
                    plot(Wykres')
                    hold on
                    
                    Plots.(genvarname(strcat('Wykresy_',num2str(ProblematyczneParyCounter))))(ProblematyczneParyCounter,:,:)=Wykres;% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU
                                                      
                    figure (3)
                    hold on                    
                    plot(WykresAvg', 'Color',ColorScale(OdKtorego(MovieNumberCounter),:))
                    hold on
                    
                    ProblematyczneParyCounter=ProblematyczneParyCounter+1;
                
                 end          
                    
                Hist.(genvarname(strcat('Movie_',num2str(MovieNumberCounter))))(LicznikParElektrod)=OdKtorego(MovieNumberCounter);% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU
                
        end

    LicznikParElektrod=LicznikParElektrod+1;
    end
end

getProblematycznePary= ProblematycznePary(:,1:ProblematyczneParyCounter-1)

for HistCounter=1:NoOfHistograms     %%%%%%1:25   %%%generowanie wykresow do histogramu 
    figure (2)
    subplot(5,5,HistCounter);    
    hist( Hist.(genvarname(strcat('Movie_',num2str(HistCounter)))),[1:1:40]);    
    %hist(['Hist.Movie_',num2str(HistCounter)],[1:1:40]); 
    title(['Movie: ',num2str(HistCounter)],'FontSize', 10, 'FontName', 'Arial CE');
    axis([7 40 0 180]);    
    hold on
    subplot(5,5,21)
    xlabel('Nr probki od której wartoœæ jest poni¿ej cutoff')
    ylabel('Ilosc powtorzen')   
end

toSubplot=ceil((ProblematyczneParyCounter-1)/2);
% for ProblematyczneParyCounter=1:ProblematyczneParyCounter-1
%     figure (4)
%     hold on
%     subplot(toSubplot,2,ProblematyczneParyCounter)
%     title(['Stim. el: ',num2str(ProblematycznePary(1,ProblematyczneParyCounter)),'\newlineRead el. : ', num2str(ProblematycznePary(2,ProblematyczneParyCounter))])
%     hold on
%     %     
%     size(Plots.Wykresy_1)
%     Nazwa=strcat('Plots.Wykresy_',num2str(ProblematyczneParyCounter))
%     toPlot=genvarname(Nazwa)
%     reshape(Nazwa,[40,23])  %size(Wykres'))
%    % reshape(Plots.Wykresy_1,40,23)
%     hold on
% end

figure(3)
title(['stim. el. avg-plots \newline \newline OdKtorego > ',num2str(OdKtorego_Analisys)])






