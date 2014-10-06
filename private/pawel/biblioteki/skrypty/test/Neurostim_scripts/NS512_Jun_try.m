TimeShiftInMs=7.5; %change
DelayInMs=5; %change

TimeShift=TimeShiftInMs*20; %time shift in sampling periods (do not change)
Delay=DelayInMs*20; % in sampling periods (do not change)

NumberOfSamples=10000; % do not change
electrodes=[1:512]; %do not change
Array=eye(512,512); % do not change
Times=[];

%MOdify this block:
El=35;
%El=[StartEl:StartEl+63];
for i=1:1
    %Array(StartEl-1+i,i)=1;
    Times=[Times Delay+(i-1)*TimeShift];
end
%Do not modify thing below
Patterns=El;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns); %do not change
MovieChunksFile=[1 Chunk]; % do not change
%break;
cd H:\2012_11\Caltech2012\rozne\testowanie_systemu\stim_files; %change
fid = fopen('el35_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('el35_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('el35_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 