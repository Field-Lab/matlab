function ClusterFileName=NS_CreateClusterFile(WritePath,FileName,NumberOfMovies,NumberOfPatterns,NumberOfRepetitions);
%Creates the cluster file that will be used later for saving the clustering
%results. 

ClusterFileName=[WritePath filesep 'ClusterFile_' FileName];
fid=fopen(ClusterFileName,'wb','ieee-le');
a=ones(NumberOfMovies*NumberOfPatterns*NumberOfRepetitions+3,1);
a(1)=NumberOfMovies; %126
a(2)=NumberOfPatterns; %330
a(3)=NumberOfRepetitions;
fwrite(fid,a,'int16');
fclose(fid);