function [TracesAll Art var0 listAmps listCurrents stimElecs firstArt]=loadTracesNoSubtract(pathToAnalysisData,patternNo,Tmax,nTrial,subSampleRate,varargin)
%load data from a pattern in a given folder
%also loads stimulus data and construct a firt artifact estimate (means
%accross trials). nTrial is the maximum number of trials
%It authomatically sorts traces increasingly with amplitude
%Gonzalo Mena,3/2016
Res=NaN;
nTrial=250;
findBundle=0;

if(nargin==6)
    params=varargin{1};
    Tmax=params.global.Tmax;
    findBundle=params.bundle.findBundle;
    findBundleTimes=params.bundle.findBundleTimes;
    nNeighborsExclude=params.bundle.nNeighborsExclude;
    detectThreshold= params.bundle.detectThreshold;
    
end

movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);



TracesAll=NaN*zeros(length(movieNos),nTrial,512,Tmax);

nTrials=zeros(1,length(movieNos));

for m=1:length(movieNos);
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
        movieNos(m), 99999);
    
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
      listStimElecs(m,:)         = channelsWithStim;
    Amp=abs(amps(1));
    
    if(m==1)
        nTrials(m)=size(dataTraces,1);
        listAmps(1)=Amp;
        listCurrents(1)=currentRangesUsed(1);
        TracesAll(1,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
        
        continue
    end
    
    
    findAmp=find(Amp==listAmps);
    
    if(isempty(findAmp))
        aux=[listAmps Amp];
        [a b]=sort([listAmps Amp],'ascend');
        
        
        index=find(b==length(aux));
        if(index<length(aux))
            
            TracesAllOld=TracesAll(1:length(listAmps),:,:,:);
            
            TracesAll(index,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
            
            TracesAll(setdiff([1:length(aux)],index),:,:,:)=TracesAllOld;
            nTrialsOld=nTrials(1:length(listAmps));
            nTrials(index)=size(dataTraces,1);
            nTrials(setdiff([1:length(aux)],index))=nTrialsOld;
            
            listAmps=sort([listAmps,Amp],'ascend');
            listCurrentsOld=listCurrents;
            listCurrents(index)=currentRangesUsed(1);
            listCurrents(setdiff([1:length(aux)],index))=listCurrentsOld;
        else
            
            TracesAll(index,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
            nTrialsOld=nTrials(1:length(listAmps));
            nTrials(index)=size(dataTraces,1);
            nTrials(setdiff([1:length(aux)],index))=nTrialsOld;
            
            listAmps=sort([listAmps,Amp],'ascend');
            listCurrentsOld=listCurrents;
            listCurrents(index)=currentRangesUsed(1);
            listCurrents(setdiff([1:length(aux)],index))=listCurrentsOld;
        end
        
    else
        TracesAll(findAmp,nTrials(findAmp)+1:nTrials(findAmp)+size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
        nTrials(findAmp)=nTrials(findAmp)+size(dataTraces,1);
    end
    
end

nTrialsSub=floor(nTrials*subSampleRate);

stimElecs=unique(listStimElecs);
TracesAllOld=TracesAll(1:length(listAmps),1:max(nTrials),:,:);
TracesAll=NaN*zeros(length(listAmps),max(nTrialsSub),512,Tmax);


sampledTrials=NaN*zeros(length(listAmps),max(nTrials));
for m=1:length(listAmps)
    sample=sort(randsample(nTrials(m),nTrialsSub(m)));
     sampledTrials(m,sample)=1;
    if(m==1)
       firstArt=nanmean(TracesAllOld(1,sample,:,:),2);
    end
   
    TracesAll(m,1:nTrialsSub(m),:,:)=TracesAllOld(m,sample,:,:);
    
    
    if(m<=5)
        a=TracesAll(m,sample,setdiff([1:512],stimElecs(1)),:);
        
        varm(m)=nanvar(a(:));
    end
    Art(m,:,:)=squeeze(nanmean(TracesAll(m,:,:,:)));
    
end

var0=nanmean(varm(1:5));

if(findBundle)
    
    [Res]=ResidualsElectrodeSimple(Art,stimElecs(1),findBundleTimes);
    
    
    
    
    pats=getNeighbors(stimElecs(1),nNeighborsExclude);
    nConds=size(Art,1);
    
    for k=1:nConds
        
        [h p]=kstest((log(Res(k,setdiff([1:512],pats)))-nanmean(log(Res(k,setdiff([1:512],pats)))))./nanstd((log(Res(k,setdiff([1:512],pats))))));
        
        pval(k)=log(p);
    end
    
    pval(1)=NaN;
    
    aux=find(pval>detectThreshold);
    comp=lastConnectedComponent(aux);
    
    if(isempty(aux));
        onset=listAmps(2);
        onsetC=2;
        return
    else
        if(length(comp)==1&&comp(end)==nConds)
            if(aux(1)==2)
                onset=NaN;
                onsetC=NaN;
                return
            else
                condi=2;
            end
        elseif(length(comp)==1&&comp(end)<nConds);
            condi=comp(end)+1;
        else
            if(comp(end)==nConds)
                condi=comp(end-1)+1;
            else
                condi=comp(end)+1;
            end
        end
        
        onset=listAmps(condi);
        onsetC=condi;
    end
else
    onset=NaN;
    onsetC=NaN;
    pval=NaN;
end

