s=reshape(NeuronData(14,:,:),50,10000);
h=plot(s')
set(h,'Color','b')

FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\RawDataPieces\proba1'];
fid=fopen(FullName,'wb','ieee-le');                                    
fwrite(fid,NeuronData,'int32');
fclose(fid);

fid=fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\RawDataPieces\proba1','r','ieee-le'); 
a=fread(fid,'int32');
fclose(fid)
%l=length(a);
%b=reshape(a,4,l/4);