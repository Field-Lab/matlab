function [breakRecElecs]  = findBreakRecElecs(breakStimElecs,recElecs,listStimElecs)
%findBreakElecs finds breakpoints in the recording electrodes, matching the information in stimulating electrodes and recording electrodes
%               It is made the implicit assumption that breakpoints are a phenomenon observed only in stimulating electrodes, and not their neighbors.
% Input:   -breakStimElecs. A cell array (with dimensions equal to the number of stimulating electrodes) such that breakStimElecs{e} 
%                           is a vector (possible empty) containing breakpoint information for the stimulating electrodes
%          -recElecs:       vector with recording electrodes
%          -listStimElecs:  vector with estimulating electrodes
% Output:  breakRecElecs: breakpoints in the recording electrodes, with the same format as breakStimElecs
% Gonzalo Mena 6/2015 
for e=1:length(recElecs)
    
    matchElec=find(recElecs(e)==listStimElecs(1,:));
    
    if(isempty(matchElec))
        breakRecElecs{e}=[];
    else
        breakRecElecs{e}=breakStimElecs{matchElec}';
    end
end
        
