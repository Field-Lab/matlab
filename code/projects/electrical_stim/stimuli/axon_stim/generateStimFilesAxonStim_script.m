% two random orders in 2 movie chunks?

% Stimulates at the center electrode and all nearest neighbor pairs (6) as
% well as all pairs of nearest neighbors (6)

% Delays the start of the first stimulus.  Almost always set to 0.
TimeShiftInMs=0;

% Center electrodes are picked to stimulate axons.  Pick low threshold
% electrodes ideally with a neighboring electrode also stimulating the
% axon. Good distance from initial segment.

CenterElectrodes=[1 14]; %limited to 2 centers
NumberOfClusters=length(CenterElectrodes);


DelayInMs = 30; %set trigger to 2 s

NumberOfSamples=40000;

%%
[electrodes Array] = generatePatternsAxonStim(CenterElectrodes, 1);


%%
MovieChunksFile = generateMovieAxonStim(NumberOfClusters,Array,TimeShiftInMs,DelayInMs,NumberOfSamples);

keyboard

fid = fopen('axon2_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('axon2_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('axon2_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 