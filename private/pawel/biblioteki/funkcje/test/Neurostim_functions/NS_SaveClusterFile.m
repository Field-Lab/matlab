function header=NS_SaveClusterFile(FilePath,PatternNumber,MovieNumber,WaveformTypes);

%1. Load the file
fid=fopen(FilePath,'r','ieee-le');
b=fread(fid,'int16');
%size(b)
fclose(fid);
header=b(1:3)
data=b(4:length(b));
ClusterIndexes=reshape(data,header(1),header(2),header(3));

%2. Write new values
ClusterIndexes(MovieNumber,PatternNumber,1:length(WaveformTypes))=WaveformTypes;
size(ClusterIndexes)
size(reshape(ClusterIndexes,header(1)*header(2)*header(3),1));
a=[header' reshape(ClusterIndexes,header(1)*header(2)*header(3),1)']';
%figure(200)
%plot(a,'bd-')
%size(a)

%3. Save the file
FilePath;
fid1=fopen(FilePath,'wb','ieee-le');
fwrite(fid1,a,'int16');
fclose(fid1);