function breakpoints=findBreakpoints(listCurrentRangesUsed)
%findBreakpoints findBreakpoints finds hardware breakpoints some electrodes
%   Input: listCurrentRangesUsed is a matrix with with rows equal to the number of movies and
%          columns equal to the number of electrodes. It comes from the currentRangesUsed output of the getStimAmps function.
%   Output: breakpoints: an array with same dimensions as the number of electrodes, each of them  (possibly empty) is a vector of
%           breakpoint conditions, defined as the last one before a new hardware regime occurs
% Gonzalo Mena 06/15
for i=1:size(listCurrentRangesUsed,2)
    
    difAmps=diff(listCurrentRangesUsed(:,i));
    breakpoints{i}=find(~(difAmps==0));
end
