function input = addAxonBreakPoint(input,breakPoint)
% addAxonBreakPoint adds manually an Axonal Bundle Breakpoint to the input structure, defined as the
% last condition j for which there is no observed traveling wave following stimulation.
%  inputs:   input (structure) and breakPoint j (should be not equal to 1
%  or J)
%            
% output:  input (structure) with the 
% usage:  Use it after defining input.tracesInfo.breakRecElecs
%
% Gonzalo Mena 6/2015 


E = input.tracesInfo.E;
J = input.tracesInfo.J;
if( breakPoint == J || breakPoint == 1)
    return
else
for e=1:E
    input.tracesInfo.breakAxon{e} = breakPoint;
    input.tracesInfo.breakPoints{e}= unique(sort([input.tracesInfo.breakAxon{e} input.tracesInfo.breakRecElecs{e}],'ascend'));
end
end

