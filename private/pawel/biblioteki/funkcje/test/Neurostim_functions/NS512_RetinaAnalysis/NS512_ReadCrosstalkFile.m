function CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,EventNumber);
fid2=fopen(CrosstalkFilePath,'r');
a=fread(fid2,'integer*2'); %dane 1-8: pierwszy przypadek 
fclose(fid2);
CrosstalkInfoFull=reshape(a,length(a)/8,8);
%CrosstalkInfoFull(:,1);
CrosstalkInfo=CrosstalkInfoFull(EventNumber,:);