function [ SigmaAvg,SigmaArray] = sigma(  MovieNumber,CenterChannel,ChannelRead,DataPath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SigmaArray(1:23)=0;
[Avg_Plots, All_Plots] = getTTX_PlotData( MovieNumber,CenterChannel,ChannelRead,DataPath);

for DT_Count=1:23
    eval(['sgm_',num2str(DT_Count),'=All_Plots(:,DT_Count)-(Avg_Plots)''',';'])
    Nazwa=strcat('sgm_',num2str(DT_Count));
    %plot(eval(Nazwa),'LineWidth',1)
    SigmaArray(DT_Count)=std(eval(Nazwa));
    %SigmaAvg=mean(SigmaArray) %% œrednia arytmetyczna
    SigmaAvg=sqrt((sum(SigmaArray.^2))/max(size(SigmaArray)));    %% œednia kwadratowa
    
end

end

