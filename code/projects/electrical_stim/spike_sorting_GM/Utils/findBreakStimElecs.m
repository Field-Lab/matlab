function breakpoints=findBreakpoints(listCurrentRangesUsed)


for i=1:size(listCurrentRangesUsed,2)
    
    difAmps=diff(listCurrentRangesUsed(:,i));
    breakpoints{i}=find(~(difAmps==0));
end
