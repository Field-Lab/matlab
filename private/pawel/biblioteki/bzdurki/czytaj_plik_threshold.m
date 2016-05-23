TresholdFilePath='I:\backup1\analysis\2012-09-27-4\stim_scan\thresholds_2';

fid1=fopen(TresholdFilePath,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
c11=reshape(b,length(b)/4,4)
NeuronInfoCoordinate=find(c11(:,1)==NeuronID)
NeuronInfo=c11(NeuronInfoCoordinate,:);