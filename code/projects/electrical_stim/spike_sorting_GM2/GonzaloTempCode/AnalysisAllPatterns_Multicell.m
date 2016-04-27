%%
preparations={'2012-09-24-3','2014-09-10-0','2014-11-05-3','2014-11-05-8','2014-11-24-2','2015-04-09-2','2015-04-14-0','2015-05-27-0'};



path{1}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3';
dire{1}={'data003','data004','data005','data006'};
path{2}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-09-10-0';
dire{2}={'data003'};
path{3}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-05-3';
dire{3}={'data003','data001','data004'};
path{4}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-05-8';
dire{4}={'data002','data003'};
path{5}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2014-11-24-2';
dire{5}={'data002'};
path{6}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-04-09-2';
dire{6}={'data002','data003'};
path{7}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-04-14-0';
dire{7}={'data001'};
path{8}='/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/2015-05-27-0';
dire{8}={'data001'};

for num=1:561
    pathToAnalysisData=Output561(num).path;
    aa=find(pathToAnalysisData=='/');
    prep=pathToAnalysisData(aa(end-2)+1:aa(end-1)-1);
    group(num)=strmatch(['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' prep], path);
    pattern(num)=Output561(num).stimInfo.patternNo;
    stimElec(num)=Output561(num).stimInfo.listStimElecs(1,1);
    nstimElec(num)=size(Output561(num).stimInfo.listStimElecs,2);
    
    neuronId(num)=Output561(num).neuronInfo.neuronIds;
end
num0=561;
for num=1:29
    pathToAnalysisData=Output2(num).path;
    aa=find(pathToAnalysisData=='/');
    prep=pathToAnalysisData(aa(end-2)+1:aa(end-1)-1);
    group(num+num0)=strmatch(['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' prep], path);
    pattern(num+num0)=Output2(num).stimInfo.patternNo;
    stimElec(num+num0)=Output2(num).stimInfo.listStimElecs(1,1);
    nstimElec(num+num0)=size(Output2(num).stimInfo.listStimElecs,2);
    neuronId(num+num0)=Output2(num).neuronInfo.neuronIds;
end
cont=561+29+1;
for num=1:216
    if(num==124)
        continue
    end
    aa=find(pathToAnalysisData=='/');
    pathToAnalysisData=Output216(num).path;
    prep=pathToAnalysisData(aa(end-2)+1:aa(end-1)-1);
    group(cont)=strmatch(['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' prep], path);
    pattern(cont)=Output216(num).stimInfo.patternNo;
    stimElec(cont)=Output216(num).stimInfo.listStimElecs(1,1);
    nstimElec(cont)=size(Output216(num).stimInfo.listStimElecs,2);
    neuronId(cont)=Output216(num).neuronInfo.neuronIds;
    cont=cont+1;
end



nstimElecwr=nstimElec(setdiff([1:805],reps(:,2)));

for g=1:8
    indwr{g}=find(groupwr==g);
    neuronswr{g}=unique(neuronIdwr(indwr{g}));
end

shift=[10 0 0 0 0 0 0 0];
for g=1:8
    for n=1:length(neuronswr{g})
        templates{g}{n}=templates{g}{n}(:,shift(g)+1:end);
    end
end




load([pathToAnalysisData name]);
if(length(aux)>0)
    
    if(isequal(dires(i).name(2:end),num2str(patternNo)))
        path=pathAux;
    end
end


reps=[];
for i=1:length(pattern)
    
    pat=pattern(i);
    neu=neuronId(i);
    gro=group(i);
    aux=[pat neu gro];
    for j=i+1:length(pattern)
        pat2=pattern(j);
        neu2=neuronId(j);
        gro2=group(j);
        if(isequal(aux,[pat2 neu2 gro2]))
            reps=[reps;[i j pat neu gro]];
        end
    end
end


patternwr=patternind(setdiff([1:805],reps(:,2)),:);

groupwr=group(setdiff([1:805],reps(:,2)));
indexes=[[1:561 [1:29] [1:123],[125:216]]];
for i=1:length(patternwr)
    num=patternwr(i,2);
    if(num<=561)
        pathToAnalysisData=Output561(indexes(num)).path;
        RespAlg1{i}=Output561(indexes(num));
    elseif(num>=562&&num<=561+29)
        pathToAnalysisData=Output2(indexes(num)).path;
        RespAlg1{i}=Output2(indexes(num));
    else
        pathToAnalysisData=Output216(indexes(num)).path;
        RespAlg1{i}=Output216(indexes(num));
    end
    
    aa=find(pathToAnalysisData=='/');
    pathToAnalysisData=['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' pathToAnalysisData(aa(end-2)+1:end)];
    name=['elecResp_n' num2str(neuronId(num)) '_p' num2str(patternwr(i,1)) '.mat'];
    load([pathToAnalysisData name]);
    
    LatTrue{i}=elecResp.analysis.latencies;
    probTrue{i}=elecResp.analysis.successRates;
    tam(i)=length(elecResp.analysis.latencies);
    
end
for g=1:8
    indwr{g}=find(groupwr==g);
end
neuronIdwr=neuronId(setdiff([1:805],reps(:,2)));

p=17;
%% SPIKE SORTING
patternI=[474 351 3 292 367 317 140 99];
for g=1:8
   %for g=6
    %for g=1
    %to avoid undesided 'aliasing' effects.
    
    IndFolderInitial=dire{g};  %look for the pattern in data003,data004,data005,data006. If empty, will look at all folderes
    pathToPreparationInitial=[path{g} '/'];
    
    params=InitializeArray(pathToPreparationInitial,1,patternI(g));
    neuronIds=neuronswr{g};
    if(g==6)
        %a1=intersect(find(nstimElecwr==1),find(groupwr==6));
        pat=unique(patternwr(indwr{g},1));
    else
        pat=unique(patternwr(indwr{g},1));
    end
    IndFolder=dire{g};  %look for the pattern in data003,data004,data005,data006. If empty, will look at all folderes
    %it is optional, but if not stated, it may use
    
    pathToPreparation=pathToPreparationInitial;
    %patss=unique(patsind);
    for p=1:length(pat)
        patternNo=pat(p);
        patind=find(patternNo==patternwr(:,1));
        patind=intersect(find(nstimElecwr==1),intersect(patind,find(groupwr==g)));
        %patternNo1=RespAlg1{patind(1)}.stimInfo.listStimElecs(1);
        
        %information of an undesired folder.
        
        params.global.sortData=1;
        params.global.nTrial=80;
         params.global.subSampleRate=1;
         params.bundle.findBundle=1;
        params.bundle.useBundleAlg=0;
        params.global.useStimElec=1;
        params.global.useStimElectrodeBeforeBundle=1;
        params.global.useStimElectrodeAfterBundle=0;
        tic
        %[Output]=DoSpikeSortingLargeScaleNOEi(pathToPreparation,pathToEi,patternNo,neuronIds,params,templates2{g},IndFolder);
        
        %RespStim{g}{p}=Output;
        %RespStim{g}{p}.times=toc;
        tic
        params.bundle.useBundleAlg=1;
        params.global.useStimElec=1;
        params.global.useStimElectrodeBeforeBundle=1;
        params.global.useStimElectrodeAfterBundle=0;
          [Output]=DoSpikeSortingLargeScaleNOEi(pathToPreparation,pathToEi,patternNo,neuronIds,params,templates2{g},IndFolder);
        RespStimBundle{g}{p}=Output;
        RespStimBundle{g}{p}.times=toc;
    end
%         params.bundle.useBundleAlg=1;
%         params.bundle.updateFreq=1;
%         par.global.useStimElec=1;
%         tic
%         [Output]=DoSpikeSortingLargeScaleNOEi(pathToPreparation,pathToEi,patternNo,neuronIds,params,templates2{g},IndFolder);
%         RespAllFastBundleStim{g}{p}=Output;
%         RespAllFastBundleStim{g}{p}.times=toc;
%         
%         
    end
end

%% pLOT
g=3;
SimTem=CorrTemplates(templates{g},RespAll{g}{1}.neuronInfo.ActiveElectrodesAll);
rects{1}=[5 5];
rects{2}=[3 2];
rects{3}=[5 4];
rects{4}=[5 4];
rects{5}=[3 3];
rects{6}=[4 4];
rects{7}=[3 3];
a1=find(groupwr==g);

for num=1:12
    h=figure(num);
    patternNo=RespAll{g}{num}.stimInfo.patternNo;
    b1=intersect(a1,find(patternwr(:,1)==patternNo));
    pathToAnalysisData=RespAlg1{b1(1)}.path;
    aa=find(pathToAnalysisData=='/');
    pathToAnalysisData=['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15' pathToAnalysisData(aa(end-2)+1:end)];
    
    tickAmp=[0 max(abs(RespAll{g}{num}.stimInfo.listAmps))/2 max(abs(RespAll{g}{num}.stimInfo.listAmps))+0.0001];
    onset=RespAll{g}{num}.bundle.onset;
    indMaxIter=find(RespAllFast{g}{num}.Log.Iter==params.global.maxIter*length(templates{g}));
    
    amps=abs(RespAll{g}{num}.stimInfo.listAmps);
    clear tickLabel
    if(isnan(onset)&&isempty(indMaxIter))
        ticks=tickAmp;
    elseif(isnan(onset)&&length(indMaxIter)>0)
        
        ticks=sort([tickAmp amps(indMaxIter)+0.001],'ascend');
    elseif(~isnan(onset)&&isempty(indMaxIter))
        ticks=sort([tickAmp onset+0.00001],'ascend');
        
    else
        ticks=sort([tickAmp onset+0.00001 amps(indMaxIter)+0.001],'ascend');
    end
    
    
    for k=1:length(ticks)
        if(ticks(k)==onset+0.00001)
            tickLabel{k}='O';
        elseif(length(find((amps(indMaxIter)+0.001)==ticks(k)))>0)
            tickLabel{k}='MI';
        else
            tickLabel{k}=num2str(floor(10*ticks(k))/10);
        end
    end
    
    
    
    
    for i=1:length(templates{g})
        subplot(rects{g}(1),rects{g}(2),i)
        plot(abs(RespAll{g}{num}.stimInfo.listAmps),nansum(RespAllFast{g}{num}.neuronInfo.spikes{i}'>1)./RespAll{g}{num}.stimInfo.nTrials,'linewidth',2,'color','red')
        hold on
        axis([0 max(abs(RespAll{g}{num}.stimInfo.listAmps))+0.0001 0 1])
        
        
        plot(abs(RespAll{g}{num}.stimInfo.listAmps),nansum(RespAllFastBundle{g}{num}.neuronInfo.spikes{i}'>1)./RespAll{g}{num}.stimInfo.nTrials,'linewidth',2,'color','blue')
        title(['Neuron ' num2str(RespAll{g}{num}.neuronInfo.neuronIds(i)) 'p ' num2str(RespAll{g}{num}.stimInfo.patternNo) ])
        %plot([RespAll{g}{num}.bundle.onset RespAll{g}{num}.bundle.onset],[0 1],'linewidth',3,'color','black');
        set(gca,'XTick',ticks)  % This automatically sets
        set(gca,'XTickLabel',tickLabel)
        set(gca,'fontsize',13)
        set(gca,'TickLength',[0.08 0.08])
        
    end
    
    for neu=1:length(b1)
        nneu=neuronIdwr(b1(neu));
        numneu=find(nneu==neurons{g});
        subplot(rects{g}(1),rects{g}(2),numneu)
        pathToAnalysisData=RespAlg1{b1(1)}.path;
        aa=find(pathToAnalysisData=='/');
        pathToAnalysisData=['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' pathToAnalysisData(aa(end-2)+1:end)];
        
        name=['elecResp_n' num2str(nneu) '_p' num2str(patternNo) '.mat'];
        load([pathToAnalysisData name]);
        
        
        %plot(abs(Output216(num).stimInfo.listAmps'),nansum(spikesMany{num}{i}'>1)./(Output216(num).tracesInfo.I+1),'linewidth',2,'color','red')
        hold on
        
        plot(abs(RespAll{g}{num}.stimInfo.listAmps),elecResp.analysis.successRates','--','linewidth',3,'color','black')
    end
    for i=1:length(templates{g})
        h2=subplot(rects{g}(1),rects{g}(2),i);
        
        pos=get(h2,'position');
        yminnew=pos(2)+pos(4)/2+pos(4)/32;
        h3=axes('position',[pos(1)+pos(3)/64 yminnew pos(3)/3 pos(4)/2-pos(4)/16]);
        mat=nanmax(abs(templates{g}{i})');
        mat(RespAllBundle{g}{num}.stimInfo.stimElec)=100;
        [~,EIm_view]   = ei2matrix(log(mat)');
        imagesc(EIm_view,[0 log(100)])
        box on
        axis 'off'
        axis square
        
        
    end
    h2=subplot(rects{g}(1),rects{g}(2),length(templates{g})+1);
    Art=RespAll{g}{num}.stimInfo.Art(:,:,:);
    ResEnd=squeeze(RespAll{g}{num}.stimInfo.Art(end,:,:));
    stimElec=RespAllBundle{g}{num}.stimInfo.stimElec;
    neighbors=getNeighbors(stimElec,3);
    ResEnd(stimElec,:)=NaN;
    [Res]=ResidualsElectrodeSimple(Art,stimElec,setdiff([1:Tmax],7));
    
    Res(:,neighbors)=NaN;
    Res(1,:)=NaN;
    plot(abs(RespAll{g}{num}.stimInfo.listAmps),log(Res))
    axis([0 max(abs(RespAll{g}{num}.stimInfo.listAmps))+0.0001 nanmin(nanmin(log(Res))) nanmax(nanmax(log(Res)))])
    
    pos=get(h2,'position');
    yminnew=pos(2)+pos(4)/2+pos(4)/32;
    h3=axes('position',[pos(1)+pos(3)/64 yminnew pos(3)/3 pos(4)/2-pos(4)/16]);
    varend=nanvar(ResEnd');
    [~,EIm_view]   = ei2matrix(log(varend));
    
    imagesc(EIm_view)
    axis 'off'
    axis square
    
    set(h,'position',[100   100    1100   700])
end


%% ERRORS
for g=setdiff([1:8],6)
    %for i=1:length(indwr{g})
    for i=1:length(indwr{g})
        nTrialsAll{g}(i)=0;
        nTrialsAllO{g}(i)=0;
        neuronIds=neurons{g};
        neuron=neuronIdwr(indwr{g}(i));
        nneu=find(neuron==neuronIds);
        pat=unique(patternwr(indwr{g},1));
        
        patternNo=patternwr(indwr{g}(i));
        patindex=find(patternNo==pat);
        pathToAnalysisData=RespAlg1{indwr{g}(i)}.path;
        aa=find(pathToAnalysisData=='/');
        pathToAnalysisData=['/Volumes/MAC OS/Research/EJBigData/Datasetsvisitjun15/' pathToAnalysisData(aa(end-2)+1:end)]
        
        name=['elecResp_n' num2str(neuron) '_p' num2str(patternNo) '.mat'];
        load([pathToAnalysisData name]);
        onsetC=RespAllFastBundle{g}{patindex}.bundle.onsetC;
        onset=RespAllFastBundle{g}{patindex}.bundle.onset;
        difFast{g}{i}=NaN*zeros(length(elecResp.analysis.latencies),max(RespAllFast{g}{patindex}.stimInfo.nTrials));
        difBundleFast{g}{i}=difFast{g}{i};
        if(~isnan(onsetC))
        difFastO{g}{i}=NaN*zeros(onsetC-1,max(RespAllFast{g}{patindex}.stimInfo.nTrials));
        else
         difFastO{g}{i}=NaN*zeros(length(elecResp.analysis.latencies),max(RespAllFast{g}{patindex}.stimInfo.nTrials));   
        end
        difBundleFastO{g}{i}=difFastO{g}{i};
        
        current{g}(i)=max(max(abs(RespAllFastBundle{g}{patindex}.stimInfo.listAmps)));
        currentO{g}(i)=onset;
        for j=1:length(elecResp.analysis.latencies)
            if(isempty(elecResp.analysis.latencies{j}))
                continue
            end
            nTrials=length(elecResp.analysis.latencies{j});
            difFast{g}{i}(j,1:nTrials)=double((elecResp.analysis.latencies{j}>0)')-double(RespAllFast{g}{patindex}.neuronInfo.spikes{nneu}(j,1:nTrials)>0);
            difBundleFast{g}{i}(j,1:nTrials)=double((elecResp.analysis.latencies{j}>0)')-double(RespAllFastBundle{g}{patindex}.neuronInfo.spikes{nneu}(j,1:nTrials)>0);
            
            nTrialsAll{g}(i)=nTrials+nTrialsAll{g}(i);
            if(j<=onsetC-1)
             nTrialsAllO{g}(i)=nTrials+nTrialsAllO{g}(i); 
              difFastO{g}{i}(j,1:nTrials)= difFast{g}{i}(j,1:nTrials);
            difBundleFastO{g}{i}(j,1:nTrials)=difBundleFast{g}{i}(j,1:nTrials);
            end
             
        end
    end
end
    
    


for g=setdiff([1:8],6)
    for i=1:length(indwr{g})
        
        neuronIds=neurons{g};
        neuron=neuronIdwr(indwr{g}(i));
        nneu=find(neuron==neuronIds);
        pat=unique(patternwr(indwr{g},1));
        
        errorFast{g}(i,1)= nansum(nansum(abs(difFast{g}{i})));
        errorFast{g}(i,2)= nansum(nansum(abs(difBundleFast{g}{i})));
        errorFast{g}(i,3)= nansum(nansum(abs(difFastO{g}{i})));
        errorFast{g}(i,4)= nansum(nansum(abs(difBundleFastO{g}{i})));
    end
end


g=2;
clear patError
aaa=unique(patternwr(indwr{g}(find(errorFast{g}(:,2)>errorFast{g}(:,1)))));
pat=unique(patternwr(indwr{g},1));
for i=1:length(aaa)
    patError(i)=find(aaa(i)==pat);
end
trials=0;
trials0=0;
errors=zeros(1,2);
errorsO=zeros(1,2);
for g=setdiff([1:8],6)
    errors=errors+nansum(errorFast{g}(:,1:2));
    errorsO=errorsO+nansum(errorFast{g}(:,3:4));
    trials=trials+nansum(nTrialsAll{g});
    trialsO=trialsO+nansum(nTrialsAllO{g});
end
