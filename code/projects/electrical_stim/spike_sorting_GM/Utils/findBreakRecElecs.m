function [breakRecElecs]  = findBreakRecElecs(breakStimElecs,recElecs,listStimElecs,varargin)


for e=1:length(recElecs)
    
    matchElec=find(recElecs(e)==listStimElecs(1,:));
    
    if(isempty(matchElec))
        breakRecElecs{e}=[];
    else
        breakRecElecs{e}=breakStimElecs{matchElec}';
    end
end
        
