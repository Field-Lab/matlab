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
javaaddpath 'C:\Users\Jojne\Desktop\SS\praktyki\Vision.jar'      %%% !!!!!! do zmiany przy innym komp
DataPath='E:\data023';   %%% !!!!!! do zmiany przy innym komp
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);     %define the electrode map - must be different than 1 for silicon probes
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;

LicznikParElektrod=1;
Radius=1;  % odleg³oœæ 60 mikronów (2-120,3-180)
OverTreshold(1:25)=0;
OdKtorego=OverTreshold(1:25);
colormap(jet(33))
%ColorScale=[0,0,0.5977;0,0,0.6953;0,0,0.7969;0,0,0.8945;0,0,0.9961;0,0.1106,0.9961;0,0.1992,0.9961;0,0.3008,0.9961;0,0.3984,0.9961;0,0.5,0.9961;0,0.5977,0.9961;0,0.6953,0.9961;0,0.7969,0.9961;0,0.8945,0.9961;0,0.9961,0.9961;0.1016,0.9961,0.8945;0.1992,0.9961,0.7969;0.3008,0.9961,0.6953;0.3984,0.9961,0.5977;0.5,0.9961,0.5;0.5977,0.9961,0.3984;0.3008,0.9961,0.6953;0.7969,0.9961,0.1992;0.8945,0.9961,0.1016;0.9961,0.9961,0;0.9961,0.8945,0;0.9961,0.7969,0;       0.9961,0.6953,0;0.9961,0.7969,0;0.9961,0.5977,0;0.9961,0.5,0;0.9961,0.3984,0;0.9961,0.3008,0;0.9961,0.1992,0;0.9961,0.1016,0;0.9961,0,0;0.8945,0,0;0.7969,0,0;0.5977,0,0;0.5,0,0];
ColorScale=jet(33);
%%%%% ILOSC POMIAROW
NoOfHistograms=25;   %% ilosc rysowanych histogramow
NoOfMovies=25;       %% liczba movie
ElectrodeOnFig=64;
ProblematycznePary=zeros(2,50);
ProblematyczneParyCounter=1;
WhenTracesUnderCutOff=8;   %% od ktorej probki ustawiamy gdy przeieg ponizej cutoff
Kropka=1;   %%dla wyznczania C1 i C2 ze wszystkich 23 przebiegow kropka=1, dla wyznaczania C1 i C2 ze œredniej kropka =0
AnalizowanaPara=zeros(2,400);

%%%%% ANALIZA MOVIE
MovieDoAnalizy=21;   % na wykresie coordinates wyœwietlone polaczenia w przypadku osi¹gniêcia wybranej wartoœci
OdKtorego_Analisys=10;  % od ktorego przypadku uznajemy wartosc OdKtorego za zbyt duza

%%%% COORDINATES  - rysowanie matrycy elektrod
PlotCoordinates(ElectrodeOnFig);
load('MyColormaps','mycmap2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
for ElectrodeArrayCounter=1:max(size(ArrayChannelRead))   %%60 elektrod (bez 4,9,25,31,57)    
   
    CenterChannel = ArrayChannelRead(ElectrodeArrayCounter);    
    ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
    ChannelsPlot_mod=getAdjacent_modified(ChannelsPlot);  %% eliminacja elektrod uznanych za martwe     
    ChannelsPlot_BezWlasnej=ChannelsPlot_mod(2:max(size(ChannelsPlot_mod))); %%wersja bez dzialania elektrody na elektrode    
           
    for j=1:max(size(ArrayChannelRead))-1     %%%%% pierwotnie: :max(size(ChannelsPlot))
        a=ArrayChannelRead(1:ElectrodeArrayCounter-1);
        b=ArrayChannelRead(ElectrodeArrayCounter+1:end);
        NewArrayChannelRead=cat(2,a,b);
        ChannelRead=NewArrayChannelRead(j); 
               
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
        
        OdKtorego(1:25)=0;
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
                Wykres=(Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt);
                Wykres=Wykres';
                Wykres(:,1:7)=0;   %% zerowanie pierwszych 7 probek        
                WykresAvg=mean(Wykres);
                
                MovieNumber=MovieNumber-1;
                
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
              [ SigmaAvg,SigmaArray] = sigma(  MovieNumber,CenterChannel,ChannelRead,DataPath );
              CutOff=WspolczynnikSigma*SigmaAvg;
               
                if find((abs(WykresAvg))>CutOff)
                    OverTresholdMatrix=find((abs(WykresAvg))> CutOff );                
                    OdKtorego(MovieNumberCounter)=max(OverTresholdMatrix); %% od ktorej probki wykres jest ponizej cutoff
                else                     
                    OdKtorego(MovieNumberCounter)=WhenTracesUnderCutOff;   %% jezeli find nic nie wykryje automatycznie ustawiane jest od 8 (pierwszej)probki
                end
                
                Elektroda_Stymulujaca=ArrayChannelRead(ElectrodeArrayCounter);   %% polaczenia elektrod na wykresie coordinates
                Elektroda_Odczytowa=ChannelRead;
                AnalizowanaPara(:,LicznikParElektrod)=[Elektroda_Stymulujaca,Elektroda_Odczytowa];   
               
                 if MovieNumber == MovieDoAnalizy
               
                     if OdKtorego(MovieNumberCounter) > OdKtorego_Analisys   
                      
                        ProblematycznePary(:,ProblematyczneParyCounter)=[Elektroda_Stymulujaca,Elektroda_Odczytowa];                    
                        figure (1)
                        hold on                    
                        arrow([electrodeMap.getXPosition(Elektroda_Stymulujaca), electrodeMap.getYPosition(Elektroda_Stymulujaca)] ,[electrodeMap.getXPosition(Elektroda_Odczytowa), electrodeMap.getYPosition(Elektroda_Odczytowa)]...
                            ,'EdgeColor',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'FaceColor',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'Width',2,...
                        'TipAngle', 25,'BaseAngle',60,'Length',20)                       
                       
                        colorbar('YTickLabel',{'8','11','14','17','21','24','27','31','34','37','40'})
                        set(gcf,'Colormap',mycmap2)
                        hold on  

                        eval(['Wykresy_',num2str(ProblematyczneParyCounter),'=Wykres',';']) %% dane do stworzenia wykresu
                        Plots.(genvarname(strcat('Wykresy_',num2str(ProblematyczneParyCounter))))(ProblematyczneParyCounter,:,:)=Wykres;% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU

%                         figure (3)
%                         hold on                    
%                         plot(WykresAvg', 'Color',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'LineWidth',2.5)
                        ProblematyczneParyCounter=ProblematyczneParyCounter+1;
                     end
                 end      
           Hist.(genvarname(strcat('Movie_',num2str(MovieNumberCounter))))(LicznikParElektrod)=OdKtorego(MovieNumberCounter);% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU                
        end
    LicznikParElektrod=LicznikParElektrod+1;
    end
end

getProblematycznePary= ProblematycznePary(:,1:ProblematyczneParyCounter-1);
getAnalizowanePary=AnalizowanaPara(:,1:LicznikParElektrod-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(['E:\Histogram',num2str(LicznikParElektrod-1),'_parElektrod.mat'], '-struct', 'Hist');   %% zapis histogramu do pliku
figure (2)
for HistCounter=1:NoOfHistograms    %%%%%%1:25   %%%generowanie wykresow do histogramu 
    
    subplot(5,5,HistCounter);
    NS_GlobalConstans=NS_GenerateGlobalConstants(61);
    [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, HistCounter,[1:64],NS_GlobalConstans);
    hist( Hist.(genvarname(strcat('Movie_',num2str(HistCounter)))),[1:1:40]);    
    title([num2str(Amplitudes),' [uA]'],'FontSize', 10, 'FontName', 'Arial CE');
    %axis([7 40 0 100] ); 
    xlim([7 40]) 
    hold on
    subplot(5,5,21)
    xlabel('Nr probki od której wartoœæ jest poni¿ej cutoff')
    ylabel('Ilosc powtorzen')   
end

figure(1)
FullName1=['E:\WszystkiePrzypadki_ErrorAnalisys_Matryca',num2str(LicznikParElektrod-1),'_parElektrod_od',num2str(OdKtorego_Analisys),'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
h=gcf;
set(gca,'LineWidth',3.5)     
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1);

figure(2)

FullName1=['E:\WszystkiePrzypadki_Histogram',num2str(LicznikParElektrod-1),'_parElektrod_WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1); 







