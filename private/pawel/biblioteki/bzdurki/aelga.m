fid = fopen('E:\2014_04\NeuronsPatterns2.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);
DataNew=Data;

RawDataPath='F:\analysis\slices\2013-12-12-3-PH\data001';
SD=size(Data);
for i=1:length(Data)    
    NeuronID=Data(1,i);
    [PrimaryElectrode,Amplitudes]=NS512_FindPrimaryElectrode(RawDataPath,'E:\2014_04\data001','001',NeuronID);
    DataNew(3,i)=PrimaryElectrode;
    i
end
break
fid = fopen('E:\2014_04\NeuronsPatternsPrimaryElectrodes.bin','wb');
fwrite(fid,DataNew,'int32');
fclose(fid);