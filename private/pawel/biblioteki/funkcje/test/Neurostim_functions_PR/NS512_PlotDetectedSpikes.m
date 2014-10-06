function quality=NS512_PlotDetectedSpikes(TracesWithoutArtifact,ChannelsPlot,WaveformTypes,CorrectedTraces,EI);

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
N = length(ChannelsPlot);
SCT=size(TracesWithoutArtifact);
figure(11);

artifactsIndex = find(WaveformTypes == 0);      %Wskazniki, ktore przebiegi sa artefaktami
spikeIndex = find(WaveformTypes == 1); 
    
for i = 1:length(ChannelsPlot)
    EInorm=EI(i,:)/(sqrt(sum(EI(i,:).^2)));
    sum(EInorm.^2);
    if ChannelsPlot(i) % nie rozumiem tego if (PH, 2010-06-07)
        ChannelTraces = TracesWithoutArtifact(:,ChannelsPlot(i),:); %przebiegi na elektrodzie o numerze Channel lub sasiadach
        ChannelTraces2D = reshape(ChannelTraces,SCT(1),SCT(3));
        %Plotting Traces and artifacts
        subplot(6,N,i), h= plot(ChannelTraces2D');  
        text(20,-80,num2str(ChannelsPlot(i)),'Fontsize',16);
        set(h(artifactsIndex),'Color','Black');
        set(h(spikeIndex),'Color','Red');
        axis([0 40 -100 50]);
        grid on;
        h23=gca;
        set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [-100:20:40]);
           
        %Plotting Traces, artifacts and traces deselected due to timing constraints;
        %{
        subplot(5,N,i+N), ddd = plot(ChannelTraces2D');
        set(ddd(find(TimeCorelSpikes==0)),'Color','Black');
        set(ddd(find(TimeCorelSpikes==1)),'Color','Red');
        set(ddd(find(TimeCorelSpikes==2)),'Color','Blue');
        axis([0 40 -100 50]);
        h23=gca;
        set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [-100:20:40]);
        grid on;
        %}
            
        %Plotting time correlated traces (correct spikes only)
        SCTC=size(CorrectedTraces);
        ChannelsPlot(i);
        ChannelTracesCorrected = CorrectedTraces(:,i,:); %przebiegi na elektrodzie o numerze Channel lub sasiadach
        ChannelTracesCorrected2D = reshape(ChannelTracesCorrected,SCTC(1),SCTC(3));
        subplot(6,N,i+N), f = plot(ChannelTracesCorrected2D');
        set(f,'Color','Red');
        %axis([0 40 -100 50]);
        grid on;
        h23=gca;
        set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [-100:20:40]);                        
            
        %Plotting mean spike
        subplot(6,N,i+2*N), g = plot(EI(i,:));
        %text(20,-80,strcat(sprintf('%0.4g',RMSPercentOfMeanSpike),'%'),'Fontsize',16);
        set(g,'Color','Blue');
        axis([0 40 -100 50]);
        grid on;
        h23=gca;
        set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [-100:20:40]);
                        
        C2=ChannelTracesCorrected2D;
        for j=1:SCTC(1)
            C2(j,:)=ChannelTracesCorrected2D(j,:)-EI(i,:);
            Corr1(j)=sum(ChannelTracesCorrected2D(j,:).*EInorm)/sqrt(sum(ChannelTracesCorrected2D(j,:).^2));
            Corr2(j)=sum(ChannelTracesCorrected2D(j,:).*EInorm)/sqrt(sum(EInorm.^2));
        end
        subplot(6,N,i+3*N), g = plot(C2');
        grid on;
        set(g,'Color','Red');
        axis([0 40 -50 50]);
        h23=gca;
        set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [-50:20:50]);  
        
        subplot(6,N,i+4*N), g = plot(Corr1');
        grid on;
        set(g,'Color','Red');        
        h23=gca;
        set(h23,'YLim',[0.7 1]);
        %set(h23,'XLim',[0 100]);
        %set(h23, 'XTick', [0:5:40]);
        set(h23, 'YTick', [0.7:0.05:1]);  
        text(5,0.86,num2str(mean(Corr1)));
        text(5,0.76,num2str(std(Corr1)));
        
        if i==1           
         SCorTr=size(CorrectedTraces);
            L1=SCorTr(1);
            a=find(Corr1>0.9);
            if length(a)>0.85*L1            
                q1=1;
            else            
                q1=0;
            end
        end
        
        TracesMins=min(ChannelTracesCorrected2D');
        EImin=min(EI(i,:));
        A1=max(TracesMins,EImin)./min(TracesMins,EImin);
        subplot(6,N,i+5*N), g = plot(A1);
        grid on;
        set(g,'Color','Red');        
        h23=gca;
        %set(h23,'YLim',[0.7 1]);       
        %set(h23, 'YTick', [0.7:0.05:1]);  
        text(5,0.86,num2str(mean(Corr2)));
        text(5,0.76,num2str(std(Corr2)));    
        
        if i==1
            a=find(A1>0.8);
            if length(a)>0.85*L1 && i==1
                q2=1;
            else
                q2=0;
            end
        end
    end
end
quality=q1*q2;