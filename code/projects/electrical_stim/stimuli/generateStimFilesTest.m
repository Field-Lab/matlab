NumberOfSamples = 10000;
% refresh rate (trigger interval) in labview must be the same as this number!!!
electrodes=1:64;
Array=eye(64);

Patterns = [ones(1,50) 2*ones(1,50)];
Times = [10:10:50*10 10:10:50*10]*20;

[Times iSort] = sort(Times);
Patterns = Patterns(iSort);


MovieChunks = 3; %3 chunks
Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
MovieChunks = [MovieChunks Chunk];

Patterns = ones(1,50);
Times = (10:10:50*10)*20;

Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
MovieChunks = [MovieChunks Chunk];

Patterns = 2*ones(1,50);
Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);

MovieChunksFile = [MovieChunks Chunk];

fid = fopen('test_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('test_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('test_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);