NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
% 1) Find primary electrode for each cell
NeuronIDs=[6 110 139 320 349 367 636 726 873 889 1026 1190 1283 1428 1445 1518 1732 1804 1880 1908 2017 2044 2089 2269 2299 2483 2677 2722 2900 3037 3218 3245 3425 3605 3622 3917 4055 4158 4310 4401 4578 4609 4865 4895 4907 5102 5162 5255 5389 5602 5633 5798 5884 5977 6064 6110 6263 6320 6425 6473 6503 7011 7025 7069 7249 7278 7324 7427 7461 7638];

RawDataPath='E:\2012-09-27-4\data000';
VisionOutputPath='E:\analysis\2012-09-27-4\data000';
%VisionOutputPath='D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data000';
EIlength=60;
PrimaryElectrodes=NeuronIDs;
EIs=zeros(length(NeuronIDs),EIlength);

EI_FilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\EIs';
fid1=fopen(EI_FilePath,'w+');
fwrite(fid1,[512 EIlength length(NeuronIDs) NeuronIDs],'integer*2');

for i=1:length(NeuronIDs)
    ID=NeuronIDs(i)
    [MinValues,MaxValues,EI]=NS512_FindEIamplitudes(RawDataPath,VisionOutputPath,'000',ID);    
    fwrite(fid1,EI,'integer*2');           
end
fclose(fid1);