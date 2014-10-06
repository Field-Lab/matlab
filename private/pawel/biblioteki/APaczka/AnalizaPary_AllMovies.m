clc
clear all
close all

set(0, 'defaultfigurecolor', 'w')
%%% ELEKTRODY %%%%%% 
ArrayChannelRead=[1:3,5:8,10:24,26:30,32:56,58:64];   %%recznie wygenerowana tablica odczytu

%%%%%% ERROR 
%CutOff =5;    %poziom odciecia
WspolczynnikSigma=3.0;  %% wartoœc dpisywaæ do 2 miejsc po przecinnku (konflikt w kolejnoœci alfabetycznej przy zapisie do pliku)
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

%%%%% ILOSC POMIAROW
NoOfHistograms=25;   %% ilosc rysowanych histogramow
NoOfMovies=25;       %% liczba movie
ElectrodeOnFig=64;
ProblematycznePary=zeros(2,500);
ProblematyczneParyCounter=1;
WhenTracesUnderCutOff=8;   %% od ktorej probki ustawiamy gdy przeieg ponizej cutoff
Kropka=1;   %%dla wyznczania C1 i C2 ze wszystkich 23 przebiegow kropka=1, dla wyznaczania C1 i C2 ze œredniej kropka =0
AnalizowanaPara=zeros(2,400);

%%%%% ANALIZA MOVIE
MovieDoAnalizy=23;   % na wykresie coordinates wyœwietlone polaczenia w przypadku osi¹gniêcia wybranej wartoœci
OdKtorego_Analisys=7;  % od ktorego przypadku uznajemy wartosc OdKtorego za zbyt duza

%%%% COORDINATES  - rysowanie matrycy elektrod
%PlotCoordinates(ElectrodeOnFig);
ColorScale=[0.600000023841858,0,0.200000002980232;0.600000023841858,0.600000023841858,1;0.650980412960053,0.301960796117783,1;1,0,1;0.733333349227905,0,0.490196079015732;0.466666668653488,0,0.156862750649452;0.200000002980232,0,0;0.600000023841858,0,0;1,0,0;0.600000023841858,0.600000023841858,0;0,0.200000002980232,0;0.278431385755539,0.290196090936661,0.0627451017498970;0.380392163991928,0.184313729405403,0.164705887436867;0.470588237047195,0.305882364511490,0.447058826684952;0.329411774873734,0.164705887436867,0.380392163991928;0.149019613862038,0.0627451017498970,0.290196090936661;0,0,0.200000002980232;0,0,0.466666668653488;0,0,0.733333349227905;0,0,1;0,0.666666686534882,1;0,1,0.666666686534882;0,1,0;0.800000011920929,1,0;1,0.400000005960465,0;];
%load('MyColormaps','inna25')
load('MyColormaps','mycmaps2')
 CenterChannel=61;
 ChannelRead=59;
 
 
 figure (2)
% set(gcf,'Colormap',mycmaps2)
  for MovieNumberCounter=1:NoOfMovies
      subplot(2,1,1)
      [~, AllPlots] = getTTX_PlotData( MovieNumberCounter,CenterChannel,ChannelRead,DataPath); 
      plot(AllPlots,'Color',ColorScale(MovieNumberCounter,:))
      hold on
  end
   colorbar('YTickLabel',{'5','10','15','20','25'},'Location','EastOutside')
   %%first artifact const
   MovieNumber = 1;  
        [~, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);
        FirstArtifact=AllPlots;
                 
        NS_GlobalConstans=NS_GenerateGlobalConstants(61);
        [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
        FirstMovieAmplitude=Amplitudes;
        %%second artifact const

        MovieNumber = 2;
        [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);
        SecondArtifact=AllPlots;         
          
        NS_GlobalConstans=NS_GenerateGlobalConstants(61);
        [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
        SecondMovieAmplitude=Amplitudes;
        
        [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude,Kropka);             
        
  
   for MovieNumberCounter=1:NoOfMovies    %%25 movies              
            
            MovieNumber = MovieNumberCounter;            
            [~, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);    
            Avg_Przebieg1=AllPlots;

            MovieNumber = MovieNumberCounter+1; 
           [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);
            Avg_Przebieg2=AllPlots;
            
            [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude,Kropka);
                       
           if( MovieNumberCounter < 4)
                                                    
                MovieNumber =MovieNumberCounter;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

                MovieNumber =MovieNumberCounter+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;         
            else                
aaa=MovieNumberCounter-2;
                MovieNumber=MovieNumberCounter-2;

                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                
                [~, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);
                FirstArtifact=AllPlots;
                  
                MovieNumber=MovieNumberCounter-1;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;   

                [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);
                SecondArtifact=AllPlots;
                
                [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude,Kropka);
               
                MovieNumber =MovieNumberCounter;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;  

                MovieNumber =MovieNumberCounter+1;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, MovieNumber,[1:64],NS_GlobalConstans);
                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;       
            end      
                %Wykres=(Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt);
                Wykres=Artefakt;
                Wykres=Wykres';
                Wykres(:,1:7)=0;   %% zerowanie pierwszych 7 probek        
                WykresAvg=mean(Wykres);
                subplot(2,1,2)
                plot(WykresAvg,'Color',ColorScale(MovieNumberCounter,:))
                hold on
   end
   
   figure(2)   
   colorbar('YTickLabel',{'5','10','15','20','25'},'Location','EastOutside')
  
   figure(2)
FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\AnalizaPar_stim',num2str(CenterChannel),'_read',num2str(ChannelRead),'.tif'];
h=gcf;
%set(gca,'LineWidth',3.5)     
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1);
