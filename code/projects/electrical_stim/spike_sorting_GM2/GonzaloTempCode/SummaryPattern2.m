function SummaryPattern2(Output,params,num,varargin)

h=figure(num);
set(h,'units','normalized');
set(h,'position',[0 0 0.2 0.3])
templates=Output.neuronInfo.templates;
neuronIds=Output.neuronInfo.neuronIds;
listAmps=abs(Output.stimInfo.listAmps);
stimElec=Output.stimInfo.stimElec;
Residuals=Output.stimInfo.Residuals;
pathToAnalysisData=[Output.path.pathToAnalysisData];
patternNo=Output.stimInfo.patternNo;
spikes=Output.neuronInfo.spikes;
nTrials=Output.stimInfo.nTrials;
cmap=distinguishable_colors(length(neuronIds));

if(nargin==5)

            Output2=varargin{2};
          
end
if(nargin>=4)
    
    onsetC=Output.bundle.onsetC;
    if(~isnan(onsetC)&&onsetC>1)
        listAmps=listAmps(1:onsetC-1);
        ArtEnd=squeeze(Output.stimInfo.Art(onsetC-1,:,:));
        Residuals=Residuals(1:onsetC-1,:);
        nTrials=Output.stimInfo.nTrials(1:onsetC-1);
        for i=1:length(templates)
            spikes{i}=spikes{i}(1:onsetC-1,:);
            if(nargin==5)
                
                spikes2{i}=Output2.neuronInfo.spikes{i}(1:onsetC-1,:);
            end
        end
        
        
    else
        if(nargin==5)
            spikes2=Output2.neuronInfo.spikes;
        end
    end
end  
    neighbors=getNeighbors(stimElec,3);
    ArtEnd(stimElec,:)=NaN;
    
    Residuals(:,neighbors)=NaN;
    Residuals(1,:)=NaN;
    
    nPlots=length(neuronIds);
    rects(1)=floor(sqrt(nPlots));
    rects(2)=floor(nPlots/rects(1))+1;
    tickAmp=[0 max(listAmps)/2 max(listAmps)+0.0001];
    onset=Output.bundle.onset;
    indMaxIter=find(Output.Log.Iter==params.global.maxIter*length(templates));
    if(nargin>=4)
    indMaxIter=intersect(indMaxIter,[1:1:onsetC-1]);
    end
    
    clear tickLabel
    if(isnan(onset)&&isempty(indMaxIter))
        ticks=tickAmp;
    elseif(isnan(onset)&&length(indMaxIter)>0)
        
        ticks=sort([tickAmp listAmps(indMaxIter)+0.001],'ascend');
    elseif(~isnan(onset)&&isempty(indMaxIter))
        ticks=sort([tickAmp onset+0.00001],'ascend');
        
    else
        ticks=sort([tickAmp onset+0.00001 listAmps(indMaxIter)+0.001],'ascend');
    end
    
    
    for k=1:length(ticks)
        if(ticks(k)==onset+0.00001)
            tickLabel{k}='O';
        elseif(length(find((listAmps(indMaxIter)+0.001)==ticks(k)))>0)
            tickLabel{k}='MI';
        else
            tickLabel{k}=num2str(floor(10*ticks(k))/10);
        end
    end
    
    
    
    
    for i=1:length(templates)
        [erfParams projectionComplete error] = erfFitter([listAmps(1:end); nansum(spikes{i}'>0)./nTrials; nTrials  ],2, -1);
 
        h2=subplot(3,8,i);
        h2.XColor=cmap(i,:);
        h2.YColor=cmap(i,:);
        box on
        axis square
        hold on
       plot(listAmps,nansum(spikes{i}'>1)./nTrials,'o','markeredgecolor',cmap(i,:),'markersize',3,'markerfacecolor',cmap(i,:))
       hold on 
      axis([0 max(listAmps)+0.0001 0 1])
        
        if(nargin==5)
              plot(listAmps,nansum(spikes{i}'>0)./nTrials,'o','markeredgecolor',cmap(i,:),'markersize',3,'markerfacecolor',cmap(i,:)')
        
        end
        neuronId=neuronIds(i);
        %plot(abs(Output.stimInfo.listAmps),nansum(spikes{i}'>1)./Output.stimInfo.nTrials,'linewidth',2,'color','blue')
       %title(['Neuron ' num2str(neuronId) 'p ' num2str(patternNo) ])
        %set(gca,'XTick',ticks)  % This automatically sets
        %set(gca,'XTickLabel',tickLabel)
        set(gca,'fontsize',13)
        set(gca,'TickLength',[0.08 0.08])
        
         
        name=['elecResp_n' num2str(neuronId) '_p' num2str(patternNo) '.mat'];
        
       
        try
            load([pathToAnalysisData name]);
             [erfParams projectionComplete error] = erfFitter([listAmps(1:end); elecResp.analysis.successRates'; nTrials  ],2, -1);
 
            if(nargin>=4)
                if(~isnan(onsetC))
                plot(listAmps,elecResp.analysis.successRates(1:onsetC-1)','-','markeredgecolor','black','markersize',3,'markerfacecolor','black')
               plot(projectionComplete(1,:),projectionComplete(2,:),'-','markeredgecolor','black','markersize',3,'markerfacecolor','black','linewidth',2)
                  
           h2.XColor=cmap(i,:);
        h2.YColor=cmap(i,:);
                else
                    plot(listAmps,elecResp.analysis.successRates','o','markeredgecolor','black','markersize',3,'markerfacecolor','black')
                hold on
                   plot(projectionComplete(1,:),projectionComplete(2,:),'-','markeredgecolor','black','markersize',3,'markerfacecolor','black','linewidth',2)
                 
           h2.XColor=cmap(i,:);
        h2.YColor=cmap(i,:);
                end
            else
                 plot(listAmps,elecResp.analysis.successRates','o','markeredgecolor','black','markersize',3,'markerfacecolor','black')
              hold on
               plot(projectionComplete(1,:),projectionComplete(2,:),'-','markeredgecolor','black','markersize',3,'markerfacecolor','black','linewidth',2)
                 
           h2.XColor=cmap(i,:);
        h2.YColor=cmap(i,:);
            end
            
           plot(projectionComplete(1,:),projectionComplete(2,:),'-','color',cmap(i,:),'linewidth',2)
        
       
          
        end
         h2.XColor=[0 0 0];
        h2.YColor=[0 0 0];
        ax2=h2;
       % ax1_pos = h2.Position; % position of first axes
set(ax2,'Position',h2.Position,...
    'XColor',cmap(i,:),'YColor',cmap(i,:),'color',[1 1 1],'linewidth',2);
box on
axis square
h2.XTick=[0 2 4.39];
h2.YTick=[0 0.5 1];
xlabel('stimulus (\muA)')
ylabel('Probability')
    end
        
        
       
        set(h,'position',[0 0 0.9 0.4])
end
    
