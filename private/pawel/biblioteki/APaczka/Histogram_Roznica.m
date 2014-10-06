clc
clear all
close all

title1='Histogram3422_parElektrod.mat';   %%% zmienic nazwe wczytywanego pliku
title2='Histogram292_parElektrodWspSigma_3.mat';    %%% zmienic nazwe wczytywanego pliku

t1=sscanf(title1,'Histogram%d',[2,inf]);
t2=sscanf(title2,'Histogram%d',[2,inf]);

%DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data023'; 
DataPath='G:\data023'; 
NoOfHistograms=25;
CenterChannel=20;
for HistCounter=1:NoOfHistograms    %%%%%%1:25   %%%generowanie wykresow do histogramu 
    
    load(title1,['Movie_',num2str(HistCounter)])                      
    subplot(5,5,HistCounter);
    NS_GlobalConstans=NS_GenerateGlobalConstants(61);
    [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, HistCounter,[1:64],NS_GlobalConstans);
    hist( eval(strcat('Movie_',num2str(HistCounter))),[1:1:40]);    
    title([num2str(Amplitudes),' [uA]'],'FontSize', 10, 'FontName', 'Arial CE');    
    xlim([8 40]);
    hold on    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','g','EdgeColor','k')  
    %h12=gca;
    %set(h12,'YScale','log')
end
hold on
for HistCounter=1:NoOfHistograms    %%%%%%1:25   %%%generowanie wykresow do histogramu 
    
    load(title2,['Movie_',num2str(HistCounter)])
    subplot(5,5,HistCounter);
    NS_GlobalConstans=NS_GenerateGlobalConstants(61);
    [~,Amplitudes]=NS_StimulatedChannels(DataPath,CenterChannel, HistCounter,[1:64],NS_GlobalConstans);
    hist( eval(strcat('Movie_',num2str(HistCounter))),[1:1:40]);    
    title([num2str(Amplitudes),' [uA]'],'FontSize', 10, 'FontName', 'Arial CE');
    xlim([8 40]);
    hold on 
 end
 subplot(5,5,21)
 xlabel('Nr probki od której wartoœæ jest poni¿ej cutoff')
 ylabel('Ilosc powtorzen')   

%zapis do pliku
%FullName1=['C:\Users\Anna\Desktop\IN¯YNIER\matlab\PRZEBIEGI_MATLAB\Roznica_Histogram_',num2str(t1),'_i_',num2str(t2),'_parElektrod.tif'];
break
FullName1=['G:\Roznica_Histogram_',num2str(t1),'_i_',num2str(t2),'_parElektrod.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]);
print(h, '-dtiff', '-r120', FullName1); 




