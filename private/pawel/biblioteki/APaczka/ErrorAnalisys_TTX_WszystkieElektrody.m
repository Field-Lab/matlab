clc
clear all
close all

set(0, 'defaultfigurecolor', 'w')
%%% ELEKTRODY %%%%%% 
ArrayChannelRead=[1:3,5:8,10:24,26:56,58:64];

%%%%%% ERROR 
CutOff =4;    %poziom odciecia

%%%% PATHS %%%%%
javaaddpath 'C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\Vision.jar'
DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023';
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);     %define the electrode map - must be different than 1 for silicon probes
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;

LicznikParElektrod=1;
Radius=1;  % odleg³oœæ 60 mikronów (2-120,3-180)
OverTreshold(1:25)=0;
OdKtorego=OverTreshold(1:25);
%%%% ILOŒÆ POMIARÓW
NoOfMovies=25;
ElektrodyNaMatrycy=63;
NoOfHistograms=25;

%%%%% ANALIZA MOVIE
MovieDoAnalizy=11;   % na wykresie coordinates wyœwietlone polaczenia w przypadku osi¹gniêcia wybranej wartoœci       

%%%% COORDINATES
figure (1)
                
   for i=1:ElektrodyNaMatrycy
      CoordinatesX=electrodeMap.getXPosition(i);
      CoordinatesY=electrodeMap.getYPosition(i);
      plot(CoordinatesX,CoordinatesY,'.');
      hold on
   end
      axis off
      hold on

for ElectrodeArrayCounter=1:1   %max(size(ArrayChannelRead))   %%60 elektrod (bez 4,9,25,57)    
   
    CenterChannel = ArrayChannelRead(ElectrodeArrayCounter);
    ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
    
    %%%%% wersja bez rezultatu elektroda na sama siebie
        ChannelsPlot_BezWlasnej=ChannelsPlot(2:max(size(ChannelsPlot)))      
     
    %%%%%%%%%
  
     
    
for j=1:max(size(ChannelsPlot_BezWlasnej))     %%%%% pierwotnie: :max(size(ChannelsPlot))

ChannelRead=ChannelsPlot(j); 
        
for i=1:NoOfMovies    %%25 movies  
    
MovieNumber = i;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
SDT=size(DataTraces);
DataTraces=DataTraces(1:23,:,1:40);
ClusterWithArtifacts=[1:SDT(1)]';

ArtifactsMatrix=size(ClusterWithArtifacts);
Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');

Avg_Przebieg1=mean(Przebieg1');

MovieNumber = i+1; 

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
SDT=size(DataTraces);
DataTraces=DataTraces(1:23,:,1:40);
ClusterWithArtifacts=[1:SDT(1)]';

ArtifactsMatrix=size(ClusterWithArtifacts);
Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');

Avg_Przebieg2=mean(Przebieg2');

       if( i < 4)    

                MovieNumber = 1;                

                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';             
                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
                Avg_Przebieg1_fixed=mean(Przebieg1');             

                FirstArtifact=Avg_Przebieg1_fixed;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECOND ARTIFACT

                MovieNumber = 2;

              
                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,CenterChannel,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';             
                                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');
                Avg_Przebieg2_fixed=mean(Przebieg2');             

                SecondArtifact=Avg_Przebieg2_fixed;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                C2=(FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude);
                C1=FirstArtifact-((FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude))*FirstMovieAmplitude;
                
                k=i;
                
               
                MovieNumber =k;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
               [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

                MovieNumber =k+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;            

                               
        else                

                MovieNumber=i-2;

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
           

                FirstArtifact=Avg_Przebieg1_fixed;

                MovieNumber=i-1;
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

                SecondArtifact=Avg_Przebieg2_fixed;

                C2=(FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude);
                C1=FirstArtifact-((FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude))*FirstMovieAmplitude;

                MovieNumber =i;
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

  

                MovieNumber =i+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;           
                               
       end      
                Wykres=(Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt);
                Wykres(1:7)=0;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
               
                if (find((abs(Wykres))>CutOff))
                    OverTresholdMatrix=find((abs(Wykres))> CutOff );                    
                    OverTreshold(i)=OverTresholdMatrix(1);
                    SizeMatrixOverTreshold=size(OverTresholdMatrix);
                    OdKtorego(i)=OverTresholdMatrix(SizeMatrixOverTreshold(2));                        
                else
                     OverTreshold(i)=8;
                     OdKtorego(i)=8;
                end
                
                ElektrodaA=ArrayChannelRead(ElectrodeArrayCounter);   %% polaczenia elektrod na wykresie coordinates
                ElektrodaB=ChannelRead;

                if MovieNumber == MovieDoAnalizy;   
                    figure (1)
                    hold on
                    X_Coords=[electrodeMap.getXPosition(ElektrodaA),electrodeMap.getXPosition(ElektrodaB)]
                    Y_Coords=[electrodeMap.getYPosition(ElektrodaA),electrodeMap.getYPosition(ElektrodaB)]
                    hold on
                    
                    
                    switch OdKtorego(i)
                        case 8
                            plot(X_Coords,Y_Coords,'r')
                            text(electrodeMap.getXPosition(ElektrodaB),electrodeMap.getYPosition(ElektrodaB), '+','FontSize',15,'VerticalAlignment','middle','HorizontalAlignment','center','Color','r')        
                        case 30
                            plot(X_Coords,Y_Coords,'g')
                            text(electrodeMap.getXPosition(ElektrodaB),electrodeMap.getYPosition(ElektrodaB), '+','FontSize',15,'VerticalAlignment','middle','HorizontalAlignment','center','Color','g') 
                        case 40                            
                            plot(X_Coords,Y_Coords,'b') 
                             text(electrodeMap.getXPosition(ElektrodaB),electrodeMap.getYPosition(ElektrodaB), '+','FontSize',15,'VerticalAlignment','middle','HorizontalAlignment','center','Color','b') 
                        otherwise
                            plot(X_Coords,Y_Coords,'g')
                    end
                    hold on  
                    
                    figure (3)
                    hold on
                     switch OdKtorego(i)
                         case 8
                            plot(Wykres,'r')
                         case 40
                            plot(Wykres,'b')
                         otherwise
                            plot(Wykres,'Color','g')
                     end                            
                    hold on
                
                end          
                    
                Hist_.(genvarname(strcat('Movie_',num2str(i))))(LicznikParElektrod)=OdKtorego(i);% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
end

LicznikParElektrod=LicznikParElektrod+1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%generowanie wykresow do histogramu


for HistCounter=1:NoOfHistograms      %%%%%%1:25    
                                                                    
    figure (2)
    subplot(5,5,HistCounter);    
    hist( Hist_.(genvarname(strcat('Movie_',num2str(HistCounter)))),[1:1:40]);    
    title(['Movie: ',num2str(HistCounter)],'FontSize', 10, 'FontName', 'Arial CE');
    axis([7 40 0 180]);    
    hold on
    subplot(5,5,21)
    xlabel('Nr probki od której wartoœæ jest poni¿ej cutoff')
    ylabel('Ilosc powtorzen')   

end








