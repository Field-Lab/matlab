TimeShiftInMs=3.75;
DelayInMs=5;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 
%Array=zeros(512,64); 
Times=[];

StartEl=[129 385];
El=[StartEl(1):StartEl(1)+63 StartEl(2):StartEl(2)+63];
for i=1:128
    %Array(StartEl-1+i,i)=1;
    Times=[Times Delay+(i-1)*TimeShift];
end
Patterns=El;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie
%break;
cd C:\home\Pawel\nauka\Caltech\testowanie_systemu\stim_files; 
fid = fopen('TwoChips_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('TwoChips_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('TwoChips_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 