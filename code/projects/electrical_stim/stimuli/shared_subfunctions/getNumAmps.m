function numAmps = getNumAmps(startVal,endVal)
val = startVal; 
numAmps = 1; 
while val <= endVal
    val = val*1.1;
    numAmps=1+numAmps; 
end
disp(['number of amplitudes tested is ' num2str(numAmps)]); 