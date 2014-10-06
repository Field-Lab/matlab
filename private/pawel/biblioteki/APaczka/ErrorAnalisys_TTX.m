clc
clear all
close all

%%% ELEKTRODY %%%%%%
PatternNumber=58;  %%58                         % odleg³oœæ 60 mikronów
ArrayChannelRead=[27 24 30 28 31 32];

% PatternNumber=28;                           % odleg³oœæ 120 mikronów
% ArrayChannelRead=[23 21 17 27 15 30 41 31 34 35 38 36 ];
% 
% PatternNumber=15;                           % odleg³oœæ 180 mikronów
% ArrayChannelRead=[19 16 13 10 20 7 23 3 27 64 29 60 32 54 35 38 44 49];
% 

%%%%%% ERROR 
CutOff =5;
MaxError=[1:25];
OverTreshold=[1:25];
OverTreshold(1:25)=0;
IleWykraczaPozaZakresDlaWybrenegoMovie=OverTreshold(1:25)
OdKtorego=OverTreshold(1:25);
%ss=IleWykraczaPozaZakresDlaWybrenegoMovie

%%%% PATHS %%%%%
DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023';


ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;


for j=1:1           % j = 6  dla 60 mikronow
                    %     12 dla 120 mikronow
                    %     18 dla 180 mikronow
ChannelRead=ArrayChannelRead(j);
    
for i=1:25
    
MovieNumber = i;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
SDT=size(DataTraces);
DataTraces=DataTraces(1:23,:,1:40);
ClusterWithArtifacts=[1:SDT(1)]';

ArtifactsMatrix=size(ClusterWithArtifacts);
Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');

Avg_Przebieg1=mean(Przebieg1');

MovieNumber = i+1; 

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
SDT=size(DataTraces);
DataTraces=DataTraces(1:23,:,1:40);
ClusterWithArtifacts=[1:SDT(1)]';

ArtifactsMatrix=size(ClusterWithArtifacts);
Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');

Avg_Przebieg2=mean(Przebieg2');

       if( i < 4)    

                MovieNumber = 1;                

                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';             
                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
                Avg_Przebieg1_fixed=mean(Przebieg1');             

                FirstArtifact=Avg_Przebieg1_fixed;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECOND ARTIFACT

                MovieNumber = 2;

              
                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
                SDT=size(DataTraces);
                DataTraces=DataTraces(1:23,:,1:40);
                ClusterWithArtifacts=[1:SDT(1)]';             
                                
                ArtifactsMatrix=size(ClusterWithArtifacts);
                Przebieg2=DataTraces(ClusterWithArtifacts,ChannelRead,:);
                Przebieg2=(reshape(Przebieg2,ArtifactsMatrix(1),40)');
                Avg_Przebieg2_fixed=mean(Przebieg2');             

                SecondArtifact=Avg_Przebieg2_fixed;
                
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                C2=(FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude);
                C1=FirstArtifact-((FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude))*FirstMovieAmplitude;
                
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

                               
        else                

                MovieNumber=i-2;

                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                FirstMovieAmplitude=Amplitudes;
                
                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
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
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                SecondMovieAmplitude=Amplitudes;   


                [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
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
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);
                Ampl=Amplitudes;
                Artefakt=C1+C2*Ampl;

  

                MovieNumber =i+1;
                NS_GlobalConstans=NS_GenerateGlobalConstants(61);
                [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber, MovieNumber,[1:64],NS_GlobalConstans);

                NextAmpl=Amplitudes;
                NextArtefakt=C1+C2*NextAmpl;           
                               
       end      
                Wykres=(Avg_Przebieg2-Avg_Przebieg1)-(NextArtefakt-Artefakt);
                Wykres(1:7)=0;
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
               
%                 
                figure (2)
                plot(Wykres);
                hold on
                
                
                if (find((abs(Wykres))>CutOff))
                    OverTresholdMatrix=find((abs(Wykres))> CutOff )
                    
                    OverTreshold(i)=OverTresholdMatrix(1)
                    SizeMatrixOverTreshold=size(OverTresholdMatrix)
                    OdKtorego(i)=OverTresholdMatrix(SizeMatrixOverTreshold(2))
                    IleWykraczaPozaZakresDlaWybrenegoMovie(i)=SizeMatrixOverTreshold(2)
                   
                else
                    OverTreshold(i)=7;
                    OdKtorego(i)=7
                end
                
                
%                 if abs(max(Wykres)) > abs(min(Wykres))
%                     MaxError(i)=abs(max(Wykres));
%                 else
%                     MaxError(i)=abs(min(Wykres));
%                 end          
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                      

end

% if (j == 6)
%     PatternNumber=38;
%     ArrayChannelRead=[33 36 35 44 39 43];
%     j=1;
%     i=1;
%  end    


 figure (1)
hold on
hist_znormalizowany=hist(OdKtorego);
bar(hist_znormalizowany);
title('histogram: wszystkie movie dla pojedynczej elektrody')
xlabel('Ktora probka')
ylabel('Ile razy sie pojawila')
hold on
% hold on
% 
% figure (3)
% hold on
% bar(IleWykraczaPozaZakresDlaWybrenegoMovie,'LineWidth',1);
% title('Ile wykracza poza zakres dla danego Movie(n+1)-Movie(n)');
% hold on

figure (4)
hold on
bar(OdKtorego,'LineWidth',1);
title('Od Którego');
xlabel('Movies')
ylabel('Od której probki')
hold on






end

% figure (1)
% hold on
% title(['Pierwsza próbka przekracaj¹ca zakres +/-',num2str(CutOff),'. Odleglosc elektrod: 60 mikronow.'],'FontSize', 10, 'FontName', 'Arial CE')
% set(gca,'LineWidth',2.5,'Box','on')
% grid on 
% hold on
% 
% FullName=['C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\Error_TTX_',num2str(j*10),'mikronow.tif'];
% h=gcf;
% set(h,'PaperUnits','inches');
% set(h,'PaperSize',[16 9]);
% set(h,'PaperPosition',[0 0 16 9]);
% print(h, '-dtiff', '-r120', FullName);
% 







