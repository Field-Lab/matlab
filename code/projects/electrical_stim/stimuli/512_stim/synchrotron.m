%% Synchrotron electrical stimulation
% Script that activates electrode pairs sequentially in time (separated by
% 0, 2, 4, 8, 16, 32 samples). A central electrode is chosen, and the
% surrounding electrodes pulse before the center electrode. Currently
% testing 2 ratios of current amplitudes between the surrounding and center 
% elecs: 0.5: 1 and 1: 1
% L. Grosberg February 2015

saveFiles = 1; %Set to 1 to save stimulation files, 0 for testing
%% Electrode / cell choice
% Center electrode should be on an axon hillock, defined by the electrode
% that activates a given cell at the minimum threshold. 
% Surrounding electrodes pulse before the center electrode

% **only tested for 3 center elecs**

CenterElectrodes=[204 314 457]; %2015-03-09-5 
% 204 targets 3046
% 314 targets 4776 and 4776
% 457 targets 6721 and 6841

% CenterElectrodes=[213 304 445]; %2015-03-09-0 
% % 213 targets cell
% % 304 targets cell 4547

% CenterElectrodes=[474 502 414]; %2015-03-09-0 data012
% %474 targets 7097
% % 502 targets 7533
% % 414 targets 

% CenterElectrodes=[146 187 352]; % 2015-02-24-1
% CenterElectrodes=[181 6 425]; %2015-02-24-4 
% 181 targets cell 2706
% 6 targets neuron 68
% 425 targets 6243

nClusters          = length(CenterElectrodes);
timeShiftInMs      = 0;        % offset the first stimulus from the beginning movie chunk
arraySpacing       = 60;       % 30 or 60 micron spacing
delayinms          = 7; % Time between each synchrotron pulse
interPairDelayInMs = 7.5;      % clusterDelayInMs*numClusters < delayInMs
numberOfSamples    = 10000;    % 0.5s SET TO TRIGGER INTERVAL IN LABVIEW

%%
electrodes = []; 
for n = 1:nClusters
    electrodes = [electrodes getCluster512(CenterElectrodes(n))]; %#ok<AGROW>
end
firstAmps = [0.5 1]; % Try these relative amplitudes on the outer electrodes 
nAmpRatios = length(firstAmps); 
Array = zeros(length(electrodes),nClusters*(1+6*nAmpRatios+6*nAmpRatios)); 
for i = 1:nClusters
    Array(7*(i-1)+1 : 7*(i) , 18*(i-1)+ (7*(i-1)+1 : 7*(i))) = 1*diag(ones(7,1));
    Array(7*(i-1)+2 : 7*(i) , 18*(i-1)+ (7*(i)+1 : 7*(i)+6)) = 0.5*diag(ones(6,1));
    Array(7*(i-1)+2 : 7*(i) , 18*(i-1)+ (7*(i)+7 : 7*(i)+12)) = 1*diag(ones(6,1));
    Array(7*(i-1)+2 : 7*(i) , 18*(i-1)+ (7*(i)+13 : 7*(i)+18)) = 0.5*diag(ones(6,1));
    Array(7*(i-1)+1 , 18*(i-1)+ (7*(i)+7 : 7*(i)+18)) = 1; 
end
figure; imagesc(Array); colorbar; 
xlabel('Pattern number'); 
ylabel('electrode (see electrode variable for value)'); 

% if(clusterOverlapCheck512(clusterElectrodes))
%     error('Clusters overlap.  Aborting.')
% end
%% Generate movie files
% Sample delays
delays = [2 4 8 16 32]; 
movieChunks = length(delays);
movieOrder = [reshape([ones(size(2:13)); 2:13; ...
    26*ones(size(2:13)); 27:38;
    51*ones(size(2:13)); 52:63],1,[]) reshape([14:25 ;39:50 ;64:75],1,[])];
for d = 1:length(delays)
    time1 = 0:delayinms*20:delayinms*20*(36-1);
    time2 = (0:delayinms*20:delayinms*20*(36-1))+delays(d);
    time3 = delayinms*20*(36) : delayinms*20 : delayinms*20*(36+35);
    allTimes = sort(cat(2,time1,time2,time3));
    Chunk=NS_MovieChunkGenerationForExperiment(allTimes, numberOfSamples, movieOrder);
    movieChunks = [movieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

MovieChunksFile=movieChunks;
%% Save files
if saveFiles
    fid = fopen('synch_el','wb','ieee-le.l64');
    fwrite(fid,electrodes,'int32');
    fclose(fid);
    
    fid = fopen('synch_pt','wb','ieee-le.l64');
    fwrite(fid,Array,'double');
    fclose(fid);
    
    fid = fopen('synch_mv','wb','ieee-le.l64');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
end