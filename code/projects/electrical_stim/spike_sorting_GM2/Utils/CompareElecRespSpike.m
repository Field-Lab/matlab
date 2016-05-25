function [Measures Activations thresH thresAlg  nTrials timesFN params onsets]=CompareElecRespSpike(pathToAnalysisData,neuronId,patternNo,Output,varargin)


listAmps=Output.stimInfo.listAmps;
nneu=find(Output.neuronInfo.neuronIds==neuronId);
spikes=Output.neuronInfo.spikes{nneu};
        name=['elecResp_n' num2str(neuronId) '_p' num2str(patternNo) '.mat'];
        load([pathToAnalysisData name]);
        sp=NaN*zeros(size(spikes));
        elecStimAmps=abs(elecResp.stimInfo.stimAmps);
        nTrialsMax=size(spikes,2);
        
        for i=1:length(listAmps)
       ind=find(listAmps(i)==elecStimAmps);
       
       spaux=[];
       for j=1:length(ind)
           spaux=[spaux elecResp.analysis.latencies{ind(j)}'];
       end
       sp(i,1:length(spaux))=spaux;
       nTrials(i)=length(spaux);
        end
        
        Jeff0=length(find(nTrials>0));
 
        spikes0=spikes;
        sp0=sp;
   if(nargin==5)
       
       thresHolds=varargin{1};
   else
       thresHolds=-1000000000;
   end
       for p=1:length(thresHolds)
          [onset onsetC]=findBundleFrompValStandAlone(Output.bundle.pvals,listAmps,thresHolds(p));
          
          onsets(p,:)=[onset onsetC];
       
          
   if(~isnan(onset))
      
       Jeff=min(Jeff0,onsetC-1);
   else
       Jeff=Jeff0;
   end
   
   
   
        spikes=spikes0(1:Jeff,:);
        sp=sp0(1:Jeff,:);
        
        
        spikeCAlg=spikes(:);
        spikeCH=sp(:);
        
        indnonan=find(~isnan(spikeCH));
        spikeCAlg=double(spikeCAlg(indnonan));
        spikeCH=double(spikeCH(indnonan));
        spikeCAlg=spikeCAlg>0;
        latCH=spikeCH;
        spikeCH=spikeCH>0;
        
      Measures(p,:)=[length(intersect(find(spikeCAlg==1),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==1))),length(intersect(find(spikeCAlg==0),find(spikeCH==0))) length(intersect(find(spikeCAlg==1),find(spikeCH==0)))]; 
      
      timesFN{p}=latCH(intersect(find(spikeCAlg==0),find(spikeCH==1)));
Activations{p}=[nanmean(spikes'>0);nanmean(sp'>0)];
[erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff); nanmean(spikes'>0); nTrials(1:Jeff) ],2, -1);
thresAlg(p) = -erfParams(2)/erfParams(1);
params(p,1:2)=erfParams;
[erfParams projectionComplete error] = erfFitter([listAmps(1:Jeff); nanmean(sp(1:Jeff,:)'>0); nTrials(1:Jeff) ],2, -1);
thresH(p)= -erfParams(2)/erfParams(1); 
params(p,3:4)=erfParams;
       end
   

