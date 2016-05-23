fid = fopen('E:\2014_04\NeuronsPatterns2.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);

Neurons=unique(Data(1,:));
for i=1:length(Neurons)
    i
    WritePath='E:\analiza\slices\2014_04\ClusterFiles';
    ClusterFileName=NS_CreateClusterFile(WritePath,num2str(Neurons(i)),452,512,50);
end