TimeShiftInMs=10;
DelayInMs=100;

NumberOfSamples=10000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:512];
Array=eye(512); 

Patterns=[];
Times=[];
Electrodes=[1    71   141   211   281   351   421   491];

for i=1
    for el=1:length(Electrodes)
        Patterns=[Patterns Electrodes(el)];
        Times=[Times (i-1)*2000+el*150];
    end        
end

Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd C:\home\pawel\nauka\SantaCruz\Stim_system_testing; 
fid = fopen('Few_el_electrodes2','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('Few_el_patterns2','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('Few_el_movie2','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 