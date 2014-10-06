NumberOfSamples=20000; %the total length of the stimulation pattern; remember that you need to define the trigger interval in the data acquisition software to 1 seconf or longer
electrodes=[1:64]; %for simplicity: define 
Array=eye(64); % each pattern includes just one electrode; pattern number is always equal to the electrode number

Patterns=[20 21 22 23 24 20]; %patterns to use
Times=[100 2100 4100 6100 8100 12100]; %time points for all the pulses; the value is given in samples (20000 samples=1 second)

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns); %this generates the movie chunk based on two vectors - number of patterns and time points
MovieChunksFile=[1 Chunk]; %this generates the data to be saved in the "movie" file - 
% if n chunks, specify above by creating MovieChunksFile=[n Chunk]; %this generates the data to be saved in the "movie" file - 

%now just generate the files
% cd C:\home\Pawel\nauka\StimFiles\proby; %of course you need to re-define this
cd ~/Desktop/StimFiles
fid = fopen('Few_el_electrodes3','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('Few_el_patterns3','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('Few_el_movie3','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 