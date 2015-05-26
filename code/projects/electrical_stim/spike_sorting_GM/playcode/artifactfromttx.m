pathToAnalysisData = '/Users/gomena/Research/EJ-2014-11-05-Processed/data007/';


NeuronIds = [2418;2434;2913;2763];
patterns=[1443:1500];
patternNo =  1443;
elecs=[170 233 194 185 162 163];
Tmin=1;
Tmax=40;
for p=1:length(patterns)
 p
    patternNo=patterns(p);
recordingElecs      =    getElecsFromPatternFile([pathToAnalysisData 'pattern_files/'],patternNo);
elecs2{p}=recordingElecs;
movieNos           = findMovieNos([pathToAnalysisData],patternNo);

    

%[elecs newtemplatesp spikestsp tspikestsp]=CreateTemplatesfromelecResp([pathToAnalysisData 'Responses/'],p,neu,elecPats,ntrials);



    sampRate           = 20000;

    for m=1:size(movieNos,2)
       dataTraces=NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
            movieNos(m), 99999);
      
        [bbb ch]=getStimAmps([pathToAnalysisData], patternNo, movieNos(m));
     
     amps0{p}(m,:)=bbb';
     sel{p}(m,:)=ch;
    
    
         for el=1:length(elecs)
%            
            stims0{p}{m,el}=squeeze(dataTraces(:,elecs(el),1:T));
         end
    end
end


for p=11:26
    p
    patternNo=patterns(p);
 movieNos           = findMovieNos([pathToAnalysisData],patternNo);
  for m=1:size(movieNos,2)
       %dataTraces=NS_ReadPreprocessedData([pathToAnalysisData], '', 0, patternNo,...
            %movieNos(m), 99999);
      
        [bbb ch cc dd ee ff ]=getStimAmps([pathToAnalysisData], patternNo, movieNos(m));
        amps2{p}(m,:)=bbb';
        sel2{p}=ch;
        ampi0{p}(m,:)=ff;
    patternNo=patterns(p);
  end
end

for ind=1:10
for j=1:4
   for i=1:17
       
   tspikestall{ind}{j}(i,:)=[tspikest0{ind}{j}(2*i-1,:) tspikest0{ind}{j}(2*i,:)];
   stimsall{ind}{i,j}=[stims0{ind}{2*i-1,j};stims0{ind}{2*i,j}];    
   end

   spikestall{ind}{j}=tspikestall{ind}{j}>0;
   
end
amps{ind}=amps0{ind}(1:2:end,:);

[stimsr spikestr tspikestr]=SubSampleTimeTrialStims(stimsall{ind},spikestall{ind},tspikestall{ind},ones(1,17),1,40);
stimsall{ind}=stimsr;
tspikestall{ind}=tspikestr;
spikestall{ind}=spikestr;
[c ]=FindBreakpoints(amps{ind},[0.26 1.05]);
br{ind}=c;
end


for i=1:10
    ampi{i}=ampi{i}(1:2:end);
end

for ind=1:58
    breakpoints{ind}=findBreakpoints(ampi{ind});
end
for ind=1:58
    
   if(ind<=10)
   ind2=34;
   else
       ind2=17;
   end
    for j=1:ind2
    aux=[];
        for e=1:6
   
% if(ind<=10)
 %      aux2=[stims0{ind}{2*j-1,e};stims0{ind}{2*j,e}];
 %else
     aux2=[stims0{ind}{j,e};];
 %end
       aux=[aux aux2];
        end
        ArtE{ind}(j,:)=mean(aux(2:end,:));
    
       
    end
   
end
   
for i=1:58
    if(i<=10)
    amps{i}=amps0{i}(1:2:end,:);
    else
        amps{i}=amps0{i};
    end
for i=1:58
    [c ]=findBreakpoints(ampi{i});
br{i}=c;
end



    
    
for ind=7:20
figure(ind)

  for i=1:6
      indel=find(sel{ind}(1,:)==elecs(i));
      if(isempty(indel))
          brr=[];
          subplot(2,3,i)
  plotArtifact2([1:40],ArtE{ind}(:,(i-1)*40+1:(i)*40),brr,[])
  %title(num2str(br{ind}))
  xlabel(num2str(elecs))
      else
          brr=br{ind}{indel};
      
  subplot(2,3,i)
  plotArtifact2([1:40],ArtE{ind}(:,(i-1)*40+1:(i)*40),brr,[])
  %title(num2str(br{ind}))
  xlabel(num2str(elecs))
      end
      end
  end



Art=ArtE{2}(:,81:120);