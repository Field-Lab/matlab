fid = fopen('J:\2014_04\NeuronsPatterns2.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);
DataNew=Data;

PrimaryElectrodesFilePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\figures4test'
load PrimaryElectrodes;

NeuronFilePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001\data001.neurons';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates2.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);

%RawDataPath='F:\analysis\slices\2013-12-12-3-PH\data001';
SD=size(Data);
for i=1:length(Data)    
    NeuronID=Data(1,i)
    NeuronIndexInPrimaryNeuronsTable=find(PrimaryNeurons==NeuronID)
    PrimaryElectrodes(NeuronIndexInPrimaryNeuronsTable)
    %[PrimaryElectrode,Amplitudes]=NS512_FindPrimaryElectrode(RawDataPath,'E:\2014_04\data001','001',NeuronID);
    DataNew(3,i)=PrimaryElectrodes(NeuronIndexInPrimaryNeuronsTable);
    i;
end

fid = fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\figures4test\NeuronsPatternsPrimaryElectrodesNew.bin','wb');
fwrite(fid,DataNew,'int32');
fclose(fid);