%Script for TTX figures


white=[1 1 1];
cyan=[0, 153, 153]/256;
blue=[0 0 153]/256;
green=[0, 128, 0]/256;
darkBlue=[0, 77, 153]/256
darkRed= [153, 0, 0]/256
orange =[255, 153, 51]/256
yellow =[255 255 0]/256
pathToAnalysisData='/Users/gomena/Research/EJBigData/EJ-2014-11-05-Processed/data007/';
patternNo=1450;
movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);


clear data
clear dataaux
clear listCurrentRangesUsedAux
clear listAmpsAux
clear listStimElecsAux
    clear dataMeans
    clear listAmps
    clear listStimElecs
    clear listCurrentRangesUsed
    
for m = 1:length(movieNos)
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
        movieNos(m), 99999);
     [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
    
    
    listAmps(m,:)              = amps';
    listStimElecs(m,:)         = channelsWithStim;
    listCurrentRangesUsed(m,:) = currentRangesUsed;
    
    stimElecs = unique(listStimElecs(:))';
    nstimElecs = length(stimElecs);
    elecs =[];
    for e = stimElecs
        clusterElectrodes = getCluster512(e);
        elecs = [elecs clusterElectrodes];
       
    end
    elecs = unique(elecs);
    elecs = setdiff(elecs,stimElecs);
    elecs =[stimElecs elecs];

    for e = 1:length(elecs)
        
        data{m,e}=squeeze(dataTraces(:,elecs(e),1:40));
        
    end
   

end


    index1=[1:size(data,1)];
    cont=1;
    
    while(~isempty(index1))
        
        amps = listAmps(index1(1),:);
        
        
        indequal = strmatch(amps,listAmps)';
        
        
        for e=1:size(data,2)
            
            dataaux{cont,e} = [];
            
            for l=indequal
                
                dataaux{cont,e} = [dataaux{cont,e}; data{l,e}];
                
            end
            
        end
        
        listAmpsAux(cont,:)                =  listAmps(index1(1),:);
        listStimElecsAux(cont,:)            =  listStimElecs(index1(1),:);
        listCurrentRangesUsedAux(cont,:)   =  listCurrentRangesUsed(index1(1),:);
        
        index1=setdiff(index1,indequal);
        cont=cont+1;
        
    end

data                  = dataaux;
listCurrentRangesUsed = listCurrentRangesUsedAux;
listStimElecs          = listStimElecsAux;
listAmps               = listAmpsAux;

for e=1:length(elecs)
    for m=1:size(data,1)
        dataMeans{e}(m,:)=mean(data{m,e}(2:end,:));
    end
end

breakStimElecs          = findBreakStimElecs(listCurrentRangesUsed);
[breakRecElecs]         = findBreakRecElecs(breakStimElecs,elecs,listStimElecs);

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
   print(titlefig{e},'-depsc2')
   print(titlefig{e},'-djpeg')
end


titles={'Stimulating electrode 1','Stimulating electrode 2','Non-stimulating electrode'};
titlefig={'ArtTraceS1','ArtTraceS2','ArtTraceNS'};


for e=1:3
     figure(e+3)
    title(titles{e},'fontsize',18)
     
        lambda=linspace(0,1,19);
       
            for m=2:20
                m2=m-1;
            plot(time,data{10,e}(m,:)','color',colorsArtifact{4}(1,:)*(1-lambda(m2))+colorsArtifact{4}(2,:)*(lambda(m2)),'linewidth',1.5);
            xlabel('time (ms)','fontsize',18); 
            ylabel('recorded daqs','fontsize',18);
            title(titles{e},'fontsize',18)
     
            hold on
            grid('on')
        
            end
   cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/FiguresMethods')
   print(titlefig{e},'-depsc2')
   print(titlefig{e},'-djpeg')
end