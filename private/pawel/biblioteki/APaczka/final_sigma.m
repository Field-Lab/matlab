clc
clear all
close all
set(0, 'defaultfigurecolor', 'w')
%%% ELEKTRODY %%%%%% 
ArrayChannelRead=[1:3,5:8,10:24,26:30,32:56,58:64];   %%recznie wygenerowana tablica odczytu elektrod - bez 4,9,25,31,57

%%% PATHS & VARIABLES %%%%%
%javaaddpath 'C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\Vision.jar' 
%javaaddpath 'C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\Vision.jar' %%% !!!!!! do zmiany przy innym komp
 javaaddpath('C:\Users\Jojne\Desktop\SS\praktyki\Vision.jar')
%DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023';   %%% !!!!!! do zmiany przy innym komp
DataPath='E:\data023';
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);     %define the electrode map - must be different than 1 for silicon probes
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;
LicznikParElektrod=1;
Radius=1;  % odleg³oœæ 60 mikronów (2-120,3-180)  %% analiza dla najbli¿szych s¹siadów
OverTreshold(1:25)=0;
OdKtorego=OverTreshold(1:25);

%%% WYBOR COLORMAP
load('MyColormaps','mycmap2')
ColorScale=[0,1,0;1,1,0;0.749019622802734,0,0.749019622802734;0,0.500000000000000,0.500000000000000;1,0.694117665290833,0.392156869173050;0.847058832645416,0.160784319043159,0;0.500000000000000,0,0.500000000000000;0.500000000000000,1,0.500000000000000;1,0,0.500000000000000;1,1,0.500000000000000;0,0.500000000000000,1;0.800000011920929,0.600000023841858,0;0.850000023841858,0.600000023841858,0;0.899999976158142,0.600000023841858,0;0.949999988079071,0.600000023841858,0;1,0.600000023841858,0;1,0.520000040531158,0;1,0.440000027418137,0;1,0.360000014305115,0;1,0.280000001192093,0;1,0.200000002980232,0;1,0.159999996423721,0;1,0.120000004768372,0;1,0.0799999982118607,0;1,0.0399999991059303,0;1,0,0;0.925000011920929,0,0;0.850000023841858,0,0;0.774999976158142,0,0;0.699999988079071,0,0;0.625000000000000,0,0;0.550000011920929,0,0;0.474999994039536,0,0;0.400000005960465,0,0;];

%%% ILOSC POMIAROW
NoOfHistograms=25;   %% ilosc rysowanych histogramow
NoOfMovies=25;       %% liczba movie
ElectrodeOnFig=64;   %% ilosc rysowanych elektrod na matrycy
PlotCoordinates(ElectrodeOnFig); %% funkcja rysujaca matryce elektrod
ProblematycznePary=zeros(2,500);  
ProblematyczneParyCounter=1;
WhenTracesUnderCutOff=8;   %% od ktorej probki ustawiamy gdy przeieg ponizej cutoff
Kropka=1;   %%dla wyznczania C1 i C2 ze wszystkich 23 przebiegow kropka=1, dla wyznaczania C1 i C2 ze œredniej kropka =0
AnalizowanaPara=zeros(2,400);

%%%%%% ERROR 
%CutOff =5;    %poziom odciecia  - obecnie wyznaczany jest jako  wsp * sigma, gdzie sigma to odch stand.
WspolczynnikSigma=3.0;  %% wartoœc dpisywaæ do 2 miejsc po przecinnku (konflikt w kolejnoœci alfabetycznej przy zapisie do pliku)

%%%%% ANALIZA MOVIE
MovieDoAnalizy=11;   % na wykresie coordinates wyœwietlone polaczenia w przypadku osi¹gniêcia wybranej wartoœci
OdKtorego_Analisys=7;  % od ktorego przypadku uznajemy wartosc OdKtorego za zbyt duza

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
for ElectrodeArrayCounter=1:max(size(ArrayChannelRead))   %%60 elektrod (bez 4,9,25,31,57)    
   
    CenterChannel = ArrayChannelRead(ElectrodeArrayCounter);    
    ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';
    ChannelsPlot_mod=getAdjacent_modified(ChannelsPlot);  %% eliminacja elektrod uznanych za martwe     
    ChannelsPlot_BezWlasnej=ChannelsPlot_mod(2:max(size(ChannelsPlot_mod))); %%wersja bez dzialania elektrody na elektrode    
           
    for j=1:max(size(ChannelsPlot_BezWlasnej))     %%%%% pierwotnie: :max(size(ChannelsPlot))

        ChannelRead=ChannelsPlot_BezWlasnej(j);        
        MovieNumber = 1;  
        
        %%first artifact const
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
                        X_Coords=[electrodeMap.getXPosition(Elektroda_Stymulujaca),electrodeMap.getXPosition(Elektroda_Odczytowa)];                       
                        Y_Coords=[electrodeMap.getYPosition(Elektroda_Stymulujaca),electrodeMap.getYPosition(Elektroda_Odczytowa)];
                        hold on                   
                        

                        %plot(X_Coords,Y_Coords,'Color',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'LineWidth',2.5)
                        arrow([electrodeMap.getXPosition(Elektroda_Stymulujaca), electrodeMap.getYPosition(Elektroda_Stymulujaca)] ,[electrodeMap.getXPosition(Elektroda_Odczytowa), electrodeMap.getYPosition(Elektroda_Odczytowa)]...
                            ,'EdgeColor',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'FaceColor',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'Width',2,...
                        'TipAngle', 25,'BaseAngle',60,'Length',20)
                        %text(electrodeMap.getXPosition(Elektroda_Odczytowa),electrodeMap.getYPosition(Elektroda_Odczytowa), '+','FontSize',20,'VerticalAlignment','middle','HorizontalAlignment','center', 'Color',ColorScale(OdKtorego(MovieNumberCounter)-7,:))  %% '+' pojawia sie na el odczytowej    
                        colorbar('YTickLabel',{'8','11','14','17','21','24','27','31','34','37','40'})
                        set(gcf,'Colormap',mycmap2)
                       
                        hold on  

                        eval(['Wykresy_',num2str(ProblematyczneParyCounter),'=Wykres',';']) %% dane do stworzenia wykresu
                        Plots.(genvarname(strcat('Wykresy_',num2str(ProblematyczneParyCounter))))(ProblematyczneParyCounter,:,:)=Wykres;% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU

                        figure (3)
                        hold on                    
                        plot(WykresAvg', 'Color',ColorScale(OdKtorego(MovieNumberCounter)-7,:),'LineWidth',2.5)
                        ProblematyczneParyCounter=ProblematyczneParyCounter+1;
                     end
                 end      
           Hist.(genvarname(strcat('Movie_',num2str(MovieNumberCounter))))(LicznikParElektrod)=OdKtorego(MovieNumberCounter);% !! DYNAMINCZETWORZENIE STRUKTURY DO HISTOGRAMU                
        end
    LicznikParElektrod=LicznikParElektrod+1;
    end
end

figure(1)
hold on
PlotCoordinates(64)
getProblematycznePary= ProblematycznePary(:,1:ProblematyczneParyCounter-1)
getAnalizowanePary=AnalizowanaPara(:,1:LicznikParElektrod-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (2)
for HistCounter=1:NoOfHistograms    %%%%%%1:25   %%%generowanie wykresow do histogramu 
    
    subplot(5,5,HistCounter);
    NS_GlobalConstans=NS_GenerateGlobalConstants(61);
    [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, HistCounter,[1:64],NS_GlobalConstans);
    hist( Hist.(genvarname(strcat('Movie_',num2str(HistCounter)))),[1:1:40]);    
    title([num2str(Amplitudes),' [uA]'],'FontSize', 10, 'FontName', 'Arial CE');   
    xlim([8 40]);
    hold on
    
end
subplot(5,5,21)
xlabel('Nr probki od której wartoœæ jest poni¿ej cutoff')
ylabel('Ilosc powtorzen')   
%save(['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\Histogram',num2str(LicznikParElektrod-1),'_parElektrod','WspSigma_',num2str(WspolczynnikSigma),'.mat'], '-struct', 'Hist');   %% zapis histogramu do pliku
save(['E:\Histogram',num2str(LicznikParElektrod-1),'_parElektrod','WspSigma_',num2str(WspolczynnikSigma),'.mat'], '-struct', 'Hist');   %% zapis histogramu do pliku

toSubplot=ceil((ProblematyczneParyCounter-1)/2);

% figure (4)
% for ProblematyczneParyCounter=1:ProblematyczneParyCounter-1    
%     subplot(toSubplot,2,ProblematyczneParyCounter)
%     grid on
%     title(['Stim. el: ',num2str(ProblematycznePary(1,ProblematyczneParyCounter)),'\newlineRead el. : ', num2str(ProblematycznePary(2,ProblematyczneParyCounter))])
%     hold on
%     Nazwa=strcat('Wykresy_',num2str(ProblematyczneParyCounter));
%     plot(eval(Nazwa),'LineWidth',1)
%     hold on
% end

figure(1)
FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\ErrorAnalisys_Matryca',num2str(LicznikParElektrod-1),'_parElektrod_od',num2str(OdKtorego_Analisys)+1,'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
FullName2=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\figures\ErrorAnalisys_Matryca',num2str(LicznikParElektrod-1),'_parElektrod_od',num2str(OdKtorego_Analisys)+1,'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.fig'];
h=gcf;
set(gca,'LineWidth',3.5)     
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1);
print(h, '-dtiff', '-r120', FullName2);

% figure(2)
% %title('HISTOGRAM')
% FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\Histogram',num2str(LicznikParElektrod-1),'_parElektrod_WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
% FullName2=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\figures\Histogram',num2str(LicznikParElektrod-1),'_parElektrod_WspolczynnikSigma_',num2str(WspolczynnikSigma),'.fig'];
% h=gcf;
% set(h,'PaperUnits','inches');
% set(h,'PaperSize',[16 9]);
% set(h,'PaperPosition',[0 0 16 9]);
% print(h, '-dtiff', '-r120', FullName1); 
% print(h, '-dtiff', '-r120', FullName2);


figure(3)
title(['Od ',num2str(OdKtorego_Analisys)+1,' próbki'],'FontSize',15)
grid on
FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\ErrorAnalisys_AvgPlot_',num2str(LicznikParElektrod-1),'_parElektrod_od_',num2str(OdKtorego_Analisys)+1,'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
FullName2=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\figures\ErrorAnalisys_AvgPlot_',num2str(LicznikParElektrod-1),'_parElektrod_od_',num2str(OdKtorego_Analisys)+1,'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.fig'];
h=gcf;    
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1);   
print(h, '-dtiff', '-r120', FullName2);  

% figure(4)
% FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\ErrorAnalisys_AllPlots_',num2str(LicznikParElektrod-1),'_parElektrod_od_',num2str(OdKtorego_Analisys),'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.tif'];
% FullName2=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\figures\ErrorAnalisys_AllPlots_',num2str(LicznikParElektrod-1),'_parElektrod_od_',num2str(OdKtorego_Analisys),'probki_mv_',num2str(MovieDoAnalizy),'WspolczynnikSigma_',num2str(WspolczynnikSigma),'.fig'];
% h=gcf;
% set(h,'PaperUnits','inches');
% set(h,'PaperSize',[16 9]);
% set(h,'PaperPosition',[0 0 16 9]);
% print(h, '-dtiff', '-r120', FullName1);   
% print(h, '-dtiff', '-r120', FullName2);  





