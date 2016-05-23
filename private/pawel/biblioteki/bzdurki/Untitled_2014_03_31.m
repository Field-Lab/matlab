fid = fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\figures3\NeuronsPatterns2.bin'],'r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);
