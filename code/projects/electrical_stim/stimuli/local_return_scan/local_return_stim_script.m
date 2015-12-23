%% SET TRIGGER INTERVAL TO 2 SECONDS

cd('/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test')

TimeShiftInMs=0;
DelayInMs=15;

NumberOfSamples=40000;
electrodes=1:64;
Array=eye(64); 

ArrayWithReturns = zeros(64, 64);

for ii = [1:64]
    clusterTemp = getCluster(ii);
    if length(clusterTemp) == 7 %center electrode not on edge
        
        ArrayWithReturns(ii,ii) = 1;
        for jj = 2:7
            ArrayWithReturns(clusterTemp(jj),ii) = -1/6;
        end
        
    end
end

ArrayComplete = [Array ArrayWithReturns];

electrodeOrder = LaurenelectrodeOrderGenerator(0);

Patterns=[electrodeOrder electrodeOrder+64];

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;


Times = (0:Delay:Delay*(length(Patterns)-1)) + TimeShift;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

% cd C:\pawel\pliki\nauka\matlab\; 
%cd /Applications/MATLAB74/work/Lauren/stimuli/stimulus_files/; 

fid = fopen('1el_local_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('1el_local_patterns','wb','ieee-le.l64')
fwrite(fid,ArrayComplete,'double');
fclose(fid);

fid = fopen('1el_local_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 