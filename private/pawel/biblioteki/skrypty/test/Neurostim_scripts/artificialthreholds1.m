TresholdFilePath='/Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_scan/thresholds_2';

fid1=fopen(TresholdFilePath,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
c11=reshape(b,length(b)/4,4)
le=length(c11(:,4));
c11(:,4)=[0.5:1.5/(le-1):2]
fid1=fopen(TresholdFilePath,'w');
b=fwrite(fid1,c11,'double');
fclose(fid1);