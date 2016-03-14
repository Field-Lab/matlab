function [TracesAll Art var0 listAmps listCurrents ]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrials)
%load data from a pattern in a given folder
%also loads stimulus data and construct a firt artifact estimate (means
%accross trials). nTrials is the maximum number of trials


movieNos  = findMovieNos(pathToAnalysisData,patternNo);
movieNos  = sort(movieNos);



TracesAll=NaN*zeros(length(movieNos),nTrials,512,Tmax);


for m=1:length(movieNos);
    
    dataTraces = NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
        movieNos(m), 99999);
    
    [amps channelsWithStim stimAmpVectors channelsConnected elecCurrentStep currentRangesUsed] = ...
        getStimAmps(pathToAnalysisData, patternNo, movieNos(m));
    
    listAmps(m)=abs(amps(1));
    listCurrents(m,:)=currentRangesUsed(1);
    if(m==1)
        firstArt=mean(dataTraces(:,:,1:Tmax),1);
        
        
    end
    
    TracesAll(m,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax)-repmat(firstArt,size(dataTraces,1),1,1);
    TracesAll2(m,1:size(dataTraces,1),:,:)=dataTraces(:,:,1:Tmax);
    
    a=TracesAll(m,:,setdiff([1:512],patternNo),:);
    varm(m)=nanvar(a(:));
    Art(m,:,:)=squeeze(nanmean(TracesAll(m,:,:,:)));
   
   
end


var0=nanmean(varm(1:5));



    
    