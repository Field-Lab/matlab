clc
clear all
close all

%%% ELEKTRODY %%%%%%
PatternNumber=58;
%ChannelRead=47;
%ArrayChannelRead=[41 26 28 29];
%ArrayChannelRead=[41 26 28 29];
ArrayChannelRead=23

kolor=['m' 'c' 'r' 'g'];

%%%% PATHS %%%%%
DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data003';
%ClusterPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data003\ClusterFile_006';


ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;
NoOfMovies=25;

%%mozliwe ustawienia dla fig 1
  %plot(mean(Przebieg1)')
        %plot(Avg_Przebieg1)                   %% œrednia artefaktów dla wybranego ,Movies,ChannelRead,PatternNumber
        %plot(Avg_Przebieg2)                   %% œrednia artefaktów dla wybranego ,Movies+1,ChannelRead,PatternNumber
        %plot(Przebieg1)                       %% artefakty dla wybranego ,Movies,ChannelRead,PatternNumber
        %plot(Przebieg2)                       %% artefakty dla wybranego ,Movies+1,ChannelRead,PatternNumber


for ArrayChannelCounter=1:1
    ClusterPath=['C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023\ClusterFile_0',num2str(ArrayChannelRead(ArrayChannelCounter))];
    ChannelRead=ArrayChannelRead(ArrayChannelCounter);
    
    for i=1:NoOfMovies
    
        MovieNumber = i;         
        [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath);              
        Avg_Przebieg1=AvgPlots;
       
        MovieNumber = i+1; 

        [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath); 
        Avg_Przebieg2=AvgPlots;

        figure(1)
        subplot(5,5,i) 
        plot(Avg_Przebieg2-Avg_Przebieg1, kolor(ArrayChannelCounter),'LineWidth',2);      %% ró¿nica œredniej wartoœci arefatku dla Movies 1:25; ró¿nica nastêpuj¹cych po sobie Movies  
        tekst=([num2str(i+1),'-',num2str(i)])
        text(28,-5,tekst);                                 %- umieszcza tekst w punkcie (x,y) 
        axis([0, 40, -20, 20]);
        set(gca,'LineWidth',2.5,'Box','on','XTick',[],'YTick',[])
        hold on

        figure (2)

        if( i < 4)    

                MovieNumber = 1;
                
                %%first art const                
                [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath);
                FirstArtifact=AvgPlots;             
              
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                %second art const
                MovieNumber = 2;
                [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath);
                SecondArtifact=AvgPlots;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;
               
                [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude);
                
                k=i;           
                MovieNumber =k;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);

                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

                MovieNumber =k+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;

                subplot(5,5,k)

                plot((NextArtefakt-Artefakt), kolor(ArrayChannelCounter),'LineWidth',2);
                tekst=([num2str(k+1),'-',num2str(k)]);
                text(28,-5,tekst);                                 %- umieszcza tekst w punkcie (x,y) 
                axis([0, 40, -20, 20]);
                set(gca,'LineWidth',2.5,'Box','on','XTick',[],'YTick',[])              
                hold on                   
                
        else                

                MovieNumber=i-2;

                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;

%                 ClusterIndexes=NS_ReadClusterFileAll(ClusterPath); %%- wczytywanie wynikow klastrowania:  26 movies,   64 patterns   53 powtorzenia
%                 AllClusters=ClusterIndexes(MovieNumber,PatternNumber,:);
%                 %ClusterWithArtifacts=find(AllClusters==1);
                
                [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath);              
                FirstArtifact=AvgPlots;             
       
                MovieNumber=i-1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;   

%                 ClusterIndexes=NS_ReadClusterFileAll(ClusterPath); %%- wczytywanie wynikow klastrowania:  26 movies,   64 patterns   53 powtorzenia
%                 AllClusters=ClusterIndexes(MovieNumber,PatternNumber,:);
%                 %ClusterWithArtifacts=find(AllClusters==1)
                
                [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead,DataPath);              
                SecondArtifact=AvgPlots;         
                [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude)
                 
                MovieNumber =i;
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;  

                MovieNumber =i+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;

                subplot(5,5,i)
                plot(NextArtefakt-Artefakt, kolor(ArrayChannelCounter),'LineWidth',2);
                tekst=([num2str(i+1),'-',num2str(i)])
                text(28,-5,tekst);                                 %- umieszcza tekst w punkcie (x,y) 
                axis([0, 40, -20, 20]);
                set(gca,'LineWidth',2.5,'Box','on','XTick',[],'YTick',[])

                hold on

        end
                figure (3)
                subplot(5,5,i)
                plot((Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt), kolor(ArrayChannelCounter),'LineWidth',2);      %% ró¿nica œredniej wartoœci arefatku dla Movies 1:25; ró¿nica nastêpuj¹cych po sobie Movies  

                tekst=([num2str(i+1),'-',num2str(i)])
                text(28,-5,tekst);                                 %- umieszcza tekst w punkcie (x,y) 
                axis([0, 40, -20, 20]);
                set(gca,'LineWidth',2.5,'Box','on','XTick',[],'YTick',[])
                hold on             
    end
end

figure(1)
hold on
subplot(5,5,5)
legend('1','2','3','4','Location','NorthEastOutside')
hold on
subplot(5,5,21)
set(gca,'XTick',[0 10 20 30 40])
set(gca,'YTick',[-20 0 10 20])
grid on 
hold on


FullName='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\TTX_Zwykly.tif';
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName);     

figure(2)

hold on
subplot(5,5,5)
legend('1','2','3','4','Location','NorthEastOutside')
hold on
subplot(5,5,21)
set(gca,'XTick',[0 10 20 30 40])
set(gca,'YTick',[-20 0 10 20])
grid on 
hold on

FullName='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\TTX_Przewidywany.tif';
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName);     

figure(3)
hold on
subplot(5,5,5);
legend('1','2','3','4','Location','NorthEastOutside')
hold on
subplot(5,5,21)
set(gca,'XTick',[0 10 20 30 40])
set(gca,'YTick',[-20 0 10 20])
grid on 
hold on


FullName='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\TTX_Zwykly-Przewidywany.tif';
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName);








