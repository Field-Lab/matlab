function [TracesAll Art var0 listAmps listCurrents onset onsetC pval Res stimElecs]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrial,varargin)
%load data from a pattern in a given folder
%also loads stimulus data and construct a firt artifact estimate (means
%accross trials). nTrial is the maximum number of trials
%Gonzalo Mena,3/2016
findBundle=0;

if(nargin==5)
    params=varargin{1};
    Tmax=params.global.Tmax;
    nTrial=params.global.nTrial;
    findBundle=params.bundle.findBundle;
    findBundleTimes=params.bundle.findBundleTimes;
    nNeighborsExclude=params.bundle.nNeighborsExclude;
    detectThreshold= params.bundle.detectThreshold;

end

movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);



TracesAll=NaN*zeros(length(movieNos),nTrial,512,Tmax);


for m=1:length(movieNos);
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
        movieNos(m), 99999);
    
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
    
    listAmps(m)           =      abs(amps(1));
    listCurrents(m,:)     =      currentRangesUsed;
    listStimElecs(m,:)         = channelsWithStim;
    
    if(m==1)
        firstArt=mean(dataTraces(:,:,1:Tmax),1);
        
        
    end
    
    TracesAll(m,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax)-repmat(firstArt,size(dataTraces,1),1,1);
    
    
    a=TracesAll(m,:,setdiff([1:512],patternNo),:);
    if(m<=5)
    varm(m)=nanvar(a(:));
    
    
    end
   Art(m,:,:)=squeeze(nanmean(TracesAll(m,:,:,:)));
end
stimElecs=unique(listStimElecs);
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
    




    
    