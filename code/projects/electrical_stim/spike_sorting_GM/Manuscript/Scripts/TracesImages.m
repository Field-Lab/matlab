%Script for traces figures
 cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/Scripts')
   load TracesExample
time=[1:40]/20;
conditions=[9 13 14];
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
    print(['EL1' num2str(i)],'-depsc2')
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
    print('TemplateEL1','-depsc2')
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
            title(['Electrode 2, condition j=' num2str(conditions(i))],'fontsize',18)
            grid('on')
    print(['EL2' num2str(i)],'-djpeg')
    print(['EL2' num2str(i)],'-depsc2')
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
    print('TemplateEL2','-depsc2')
    end
    
   
end



