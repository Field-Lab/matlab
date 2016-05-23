fid = fopen('E:\2014_04\NeuronsPatternsPrimaryElectrodes.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);
EventData=Data(:,EventNumber);