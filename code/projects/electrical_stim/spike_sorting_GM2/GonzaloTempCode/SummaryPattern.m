function SummaryPattern(Output,params,num,varargin)

h=figure(num);
templates=Output.neuronInfo.templates;
neuronIds=Output.neuronInfo.neuronIds;
listAmps=abs(Output.stimInfo.listAmps);
stimElec=Output.stimInfo.stimElec;
Residuals=Output.stimInfo.Residuals;
pathToAnalysisData=[Output.path.pathToAnalysisData];
patternNo=Output.stimInfo.patternNo;
spikes=Output.neuronInfo.spikes;
ArtEnd=squeeze(Output.stimInfo.Art(end,:,:));
nTrials=Output.stimInfo.nTrials;
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
    
    nPlots=length(neuronIds)+1;
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
        h2=subplot(rects(1),rects(2),i);
        plot(listAmps,nansum(spikes{i}'>1)./nTrials,'linewidth',2,'color','blue')
        hold on
        axis([0 max(listAmps)+0.0001 0 1])
        
        if(nargin==5)
              plot(listAmps,nansum(spikes2{i}'>1)./nTrials,'linewidth',2,'color','red')
        end
        neuronId=neuronIds(i);
        %plot(abs(Output.stimInfo.listAmps),nansum(spikes{i}'>1)./Output.stimInfo.nTrials,'linewidth',2,'color','blue')
        title(['Neuron ' num2str(neuronId) 'p ' num2str(patternNo) ])
        set(gca,'XTick',ticks)  % This automatically sets
        set(gca,'XTickLabel',tickLabel)
        set(gca,'fontsize',13)
        set(gca,'TickLength',[0.08 0.08])
        
         
        name=['elecResp_n' num2str(neuronId) '_p' num2str(patternNo) '.mat'];
        
        
        try
            load([pathToAnalysisData name]);
            
            if(nargin>=4)
                if(~isnan(onsetC))
                plot(listAmps,elecResp.analysis.successRates(1:onsetC-1)','--','linewidth',3,'color','black')
                else
                    plot(listAmps,elecResp.analysis.successRates','--','linewidth',3,'color','black')
                end
                else
                plot(listAmps,elecResp.analysis.successRates','--','linewidth',3,'color','black')
            end
            
            
        end
        
            pos=get(h2,'position');
            yminnew=pos(2)+pos(4)/2+pos(4)/32;
            h3=axes('position',[pos(1)+pos(3)/64 yminnew pos(3)/3 pos(4)/2-pos(4)/16]);
            mat=nanmax(abs(templates{i})');
            mat(Output.stimInfo.stimElec)=100;
            [~,EIm_view]   = ei2matrix(log(mat)');
            imagesc(EIm_view,[0 log(100)])
            box on
            axis 'off'
            axis square
    end
        
        
        h2=subplot(rects(1),rects(2),nPlots);
        
        plot(listAmps,log(Residuals))
        axis([0 max(listAmps)+0.0001 nanmin(nanmin(log(Residuals))) nanmax(nanmax(log(Residuals)))])
        set(gca,'XTick',ticks)  % This automatically sets
        set(gca,'XTickLabel',tickLabel)
        set(gca,'fontsize',13)
        set(gca,'TickLength',[0.08 0.08])
        
        pos=get(h2,'position');
        yminnew=pos(2)+pos(4)/2+pos(4)/32;
        h3=axes('position',[pos(1)+pos(3)/64 yminnew pos(3)/3 pos(4)/2-pos(4)/16]);
        Varend=nanvar(ArtEnd');
        [~,EIm_view]   = ei2matrix(log(Varend));
        
        imagesc(EIm_view)
        axis 'off'
        axis square
        
        set(h,'position',[100   100    1100   700])
    end
