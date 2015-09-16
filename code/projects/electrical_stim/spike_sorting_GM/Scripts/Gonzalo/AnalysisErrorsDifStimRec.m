cd('/Users/gomena/Research/EJBigData/2015-06-algorithm-results')


path=['/Volumes/Analysis/2012-09-24-3/data003/'];
patternNo=[245,244,229,228,212];
neuronId=3427;

path=['/Volumes/Analysis/2012-09-24-3/data006/'];
patternNo=[406,435,459,490,497];
neuronId=6200;

path=['/Volumes/Analysis/2012-09-24-3/data003/'];
patternNo=[216];
neuronId=3457;

indAll=[];
indBad=[];
patAll=[];
for i=1:561
    Resp=sorted_diffStimRecElecs(i);
    if(isequal(Resp.path,path))
        if(isequal(Resp.neuronInfo.neuronIds,neuronId))
            indAll=[indAll i];
            patAll=[patAll Resp.stimInfo.patternNo];
            if(~isempty(find(patternNo==Resp.stimInfo.patternNo)))
            indBad=[indBad i];
            end
        end
    end
end

path='/Volumes/Analysis/2014-09-10-0/data003/';
patternNo=[333,342,345,350];
neuronId=5104;

indAll=[];
indBad=[];
for i=1:561
    Resp=sorted_diffStimRecElecs(i);
    if(isequal(Resp.path,path))
        if(isequal(Resp.neuronInfo.neuronIds,neuronId))
            indAll=[indAll i];
            if(~isempty(find(patternNo==Resp.stimInfo.patternNo)))
            indBad=[indBad i];
            end
        end
    end
end



neuronId=5120;
patternNo=[334,338 346];

indAll=[];
indBad=[];
for i=1:561
    Resp=sorted_diffStimRecElecs(i);
    if(isequal(Resp.path,path))
        if(isequal(Resp.neuronInfo.neuronIds,neuronId))
            indAll=[indAll i];
            if(~isempty(find(patternNo==Resp.stimInfo.patternNo)))
            indBad=[indBad i];
            end
        end
    end
end

for i=1:6
    subplot(3,2,i)
    plot(Resp.tracesInfo.data{i+k}');
end

i=7
str=['elecResp_n3457_p' num2str(patAll(i))];
load(str)
[nansum(Output(i).spikes{1}')./Output(i).tracesInfo.I;elecResp.analysis.successRates']
i=i+1;

   	for i=1:29
leg{i}=[Resp.path 'elecResp_n' num2str(Resp.neuronInfo.neuronIds) '_p' num2str(Resp.stimInfo.patternNo)]
    end

    for i=1:29
Resp=sorted_sameStimRecElecs(i);
leg{i}=[Resp.path 'elecResp_n' num2str(Resp.neuronInfo.neuronIds) '_p' num2str(Resp.stimInfo.patternNo)]
end