TimeShiftInMs=50;
DelayInMs=50;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 
%Array=zeros(512,64); 
Times=[];

El=[10:64:512];

Patterns=[];
Disconnecting=zeros(1,10);
TTLs=ones(1,5)*(-1);

Patterns=[El(1) El(2) Disconnecting El(3) TTLs El(4) Disconnecting TTLs El(5) El(6) Disconnecting El(7) TTLs El(8) Disconnecting TTLs];
TimesInMs1=[50 100 101.2:0.05:101.65 150 151.2:0.05:151.4 200 201.2:0.05:201.65 201.2:0.05:201.4];
Times=[TimesInMs1*20 (TimesInMs1+200)*20]

%Patterns=El;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie
break
cd C:\home\Pawel\nauka\StimFiles\subretinal; 
fid = fopen('Every_8th_channel_disc_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('Every_8th_channel_disc_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('Every_8th_channel_disc_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 