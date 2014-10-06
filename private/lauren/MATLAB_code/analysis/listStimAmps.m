function [allAmps, numAmps] = listStimAmps(startVal,endVal)
val = startVal;
i = 1; 
allAmps(i) = val; 
while val < endVal
    val = val*1.1; % 10% increments
    i = i+1;
    allAmps(i) = val; 
end
numAmps = length(allAmps); 
end