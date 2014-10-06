function ClusterIndex=NS_ReadClusterFile(FilePath,MovieNumber,PatternNumber,TracesNumberLimit);
%This function shows the data on the array layout. 
%FullName=[ArtifactDataPath '\' 'p' num2str(PatternNumber) '_m' num2str(MovieNumber)];

fid=fopen(FilePath,'r','ieee-le');
b=fread(fid,'int16');
fclose(fid);
header=b(1:3)
data=b(4:length(b));
%size(b)
ClusterIndexes=reshape(data,header(1),header(2),header(3));
TracesNumber=min(header(3),TracesNumberLimit);
ClusterIndex=reshape(ClusterIndexes(MovieNumber,PatternNumber,1:TracesNumber),TracesNumber,1);