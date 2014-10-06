TimeShiftInMs=7.5;
DelayInMs=5;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=zeros(512,64); 
Times=[];

el16=[1 2 3 4 9 10 11 12 17 18 19 20 25 26 27 28];
el64=[];
for i=1:16
    el64=[el64 el16(i) el16(i)+4 el16(i)+32 el16(i)+36];
end
el64=el64;
for i=1:64
    for j=0:7
        Array(el64(i)+j*64,i)=1;
    end        
    Times=[Times Delay+(i-1)*TimeShift];
end
     
Patterns=[1:64];
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie
break;
cd C:\home\pawel\praca\stim_files; 
fid = fopen('AllElectrodes_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('AllElectrodes_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('AllElectrodes_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 