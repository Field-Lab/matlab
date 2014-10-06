Frequency=60; % frequency we want to stimulate
TimeShiftInMs=1000/Frequency; %change
DelayInMs=1; %change

TimeShift=TimeShiftInMs*20; %time shift in sampling periods (do not change)
Delay=DelayInMs*20; % in sampling periods (do not change)

NumberOfSamples=10000; % do not change
electrodes=[1:512]; %do not change
Array=eye(512,512); % do not change
Times=[];

%MOdify this block:
%El=[35];
Ch=15; %no. of channel
n=30; %no. of times
El=repmat(Ch,1,n);

%El=[StartEl:StartEl+63];
for i=1:n
    %Array(StartEl-1+i,i)=1;
   Times=[Times Delay+(i-1)*TimeShift];
end
%Do not modify thing below

Patterns=El;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns); %do not change
MovieChunksFile=[1 Chunk]; % do not change
break;
cd F:\Pawel\stm_files\jun; %change
fid = fopen('el15_el_60Hz','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('el15_pt_60','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('el15_mv_60','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 