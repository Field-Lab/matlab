% 1. General constants - do not modify
NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

% 2. Define of pulse timings and widths
TTL_pulse_timings=[50:100:450]; %in miliseconds
TTL_pulse_width=1; % in miliseconds

Disconnect_timings=[149.9 349.9]; % in miliseconds
Disconnect_width=1.2; % in miliseconds

% 3. Generate the tables
Patterns=[];
Times=[];

% 3.1. Define data for TTL pulses
for i=1:length(TTL_pulse_timings)
    p1=ones(1,TTL_pulse_width*20)*(-1); % multiplication by 20 to get value in samples; -1 is the pattern number for TTL
    t1=TTL_pulse_timings(i)*20:(TTL_pulse_timings(i)+TTL_pulse_width)*20-1;    
    Patterns=[Patterns p1];
    Times=[Times t1];
end

% 3.2. Generate data for disconnection
for i=1:length(Disconnect_timings)
    p1=zeros(1,Disconnect_width*20); % multiplication by 20 to get value in samples; -1 is the pattern number for TTL
    t1=Disconnect_timings(i)*20:(Disconnect_timings(i)+Disconnect_width)*20-1;    
    Patterns=[Patterns p1];
    Times=[Times t1];
end

% 4. Generate files - do not modify
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd C:\home\Pawel\nauka\StimFiles\subretinal\stimfiles; 
fid = fopen('ex1_el','wb','l')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('ex2_pt','wb','l')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('ex3_mv','wb','l')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);