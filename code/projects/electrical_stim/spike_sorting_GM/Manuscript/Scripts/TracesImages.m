%Script for traces figures

time=[1:40]/20;
conditions=[9 13 17];
cd('/Users/gomena/Research/EJBigData/EJ-2014-11-05-Processed/data005')
load elecResp_n2418_p1459
 cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods')
   
%neuron 1
for i=1:4
    figure(i)
    if(i<=3)
        
        spikes=logical((elecResp.analysis.latencies{conditions(i)}(2:end)>0)');
        nospikes =logical(1-spikes);
        if(~all(nospikes))
            plot(time,Output.tracesInfo.data{conditions(i),1}(spikes,:)','red','linewidth',1.5)
        end
        hold on
        if(~all(spikes))
            plot(time,Output.tracesInfo.data{conditions(i),1}(nospikes,:)','blue','linewidth',1.5)
            
        end
        xlabel('time (ms)','fontsize',18);
            ylabel('recorded daqs','fontsize',18);
            title(['Electrode 1, condition j=' num2str(conditions(i))],'fontsize',18)
            grid('on')
    print(['EL1' num2str(i)],'-djpeg')
    print(['EL1' num2str(i)],'-deps')
    else
        
        plot(time,Output.neuronInfo.templates{1}([1],1:40)*1.3,'linewidth',2,'color','black')
        hold on
        plot(time,Output.neuronInfo.templates{3}([1],1:40),'--','linewidth',1.8,'color','black')
        xlabel('time (ms)','fontsize',18);
        ylabel('recorded daqs','fontsize',18);
        title('Templates recorded at electrode 1','fontsize',18)
        legend({'Neuron 1','Neuron 2'},'fontsize',16)
    grid('on')
    print('TemplateEL1','-djpeg')
    print('TemplateEL1','-deps')
    end
    
    
end


%neuron 2


cd('/Users/gomena/Research/EJBigData/EJ-2014-11-05-Processed/data005')
load elecResp_n2913_p1459
 cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods')
   

conditions=[14 15 16];
for i=1:4
    figure(i+4)
    if(i<=3)
        
        spikes=logical((elecResp.analysis.latencies{conditions(i)}(2:end)>0)');
        nospikes =logical(1-spikes);
        if(~all(nospikes))
            plot(time,Output.tracesInfo.data{conditions(i),3}(spikes,:)','red','linewidth',1.5)
        end
        hold on
        if(~all(spikes))
            plot(time,Output.tracesInfo.data{conditions(i),3}(nospikes,:)','blue','linewidth',1.5)
            
        end
        xlabel('time (ms)','fontsize',18);
            ylabel('recorded daqs','fontsize',18);
            title(['Electrode 1, condition j=' num2str(conditions(i))],'fontsize',18)
            grid('on')
    print(['EL2' num2str(i)],'-djpeg')
    print(['EL2' num2str(i)],'-deps')
    else
        
        plot(time,Output.neuronInfo.templates{3}([3],1:40)*1.3,'linewidth',2,'color','black')
        hold on
        plot(time,Output.neuronInfo.templates{1}([3],1:40),'--','linewidth',1.7,'color','black')
        xlabel('time (ms)','fontsize',18);
        ylabel('recorded daqs','fontsize',18);
        title('Templates recorded at electrode 2','fontsize',18)
        legend({'Neuron 2','Neuron 1'},'fontsize',16)
    grid('on')
     print('TemplateEL2','-djpeg')
    print('TemplateEL2','-deps')
    end
    
   
end





colors=colormap(jet);
colorsArtifact{1}=colors([1 8],:);
colorsArtifact{2}=colors([8 24],:);
colorsArtifact{3}=colors([24 40],:);
colorsArtifact{4}=colors([40 56],:);
colorsArtifact{5}=colors([56 64],:);

time=[1:40]/20;

br =breakRecElecs{1};
ranges{1}{1}=[1:br];
ranges{1}{2}=[br+1:17];
ranges{2}{1}=ranges{1}{1};
ranges{2}{2}=ranges{1}{2};
ranges{3}{1}=[1:17];
colorRanges=[2 4];
titles={'Stimulating electrode 1','Stimulating electrode 2','Non-stimulating electrode'};
titlefig={'ArtMeanS1','ArtMeanS2','ArtMeansNS'};
for e=1:3
     figure(e)
    
    for r=1:length(ranges{e})
        
        lambda=linspace(0,1,length(ranges{e}{r}));
        for m=ranges{e}{r}
            m2=m-ranges{e}{r}(1)+1;
            plot(time,dataMeans{e}(m,:),'color',colorsArtifact{colorRanges(r)}(1,:)*(1-lambda(m2))+colorsArtifact{colorRanges(r)}(2,:)*lambda(m2),'linewidth',2);
            title(titles{e},'fontsize',18)
            xlabel('time (ms)','fontsize',18); 
            ylabel('recorded daqs','fontsize',18);
            hold on
            grid('on')
        end
    end
   cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods')
   print(titlefig{e},'-deps')
   print(titlefig{e},'-djpeg')
end


titles={'Stimulating electrode 1','Stimulating electrode 2','Non-stimulating electrode'};
titlefig={'ArtTraceS1','ArtTraceS2','ArtTraceNS'};


for e=1:3
     figure(e)
    title(titles{e},'fontsize',18)
     
        lambda=linspace(0,1,19);
       
            for m=2:20
                m2=m-1;
            plot(time,data{10,e}(m,:)','color',colorsArtifact{4}(1,:)*(1-lambda(m2))+colorsArtifact{4}(2,:)*(lambda(m2)));
            xlabel('time (ms)','fontsize',18); 
            ylabel('recorded daqs','fontsize',18);
            title(titles{e},'fontsize',18)
     
            hold on
            grid('on')
        
            end
   cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods')
   print(titlefig{e},'-deps')
   print(titlefig{e},'-djpeg')
end