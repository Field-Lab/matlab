TimeShiftInMs=10;
DelayInMs=100;

NumberOfSamples=20000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:64];
Array=eye(64); 

Patterns=[];
Electrodes=[45 46];

for i=1:5
    Patterns=[Patterns Electrodes(1) Electrodes(2)];
end
        
TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;
Times=[TimeShift:Delay:TimeShift+Delay*(length(Patterns)-1)];
%Times=[TimeShift:Delay:TimeShift+Delay*60];  % good only for 61 electrodes!!

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie
break
cd C:\home\pawel\nauka\stim_files; 
fid = fopen('Few_el_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('Few_el_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('Few_el_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 