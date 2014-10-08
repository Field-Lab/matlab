%function stimulusEmulator(electrodesPath, patternsPath, moviePath, updatePeriod)

%reads stimulus files and emulates the stimulus as a movie

%% for testing
clear all

%moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/spatial_short_movie';
%patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/spatial_short_patterns';
%electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/spatial_short_electrodes';

% moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/axon_movie';
% patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/axon_patterns';
% electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/axon_electrodes';

%moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/test_movie';
%patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/test_patterns';
%electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/test_electrodes';

%moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_movie';
%patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_patterns';
%electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_electrodes';

% moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/1el_prel_movie';
% patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/1el_prel_patterns';
% electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_test/1el_prel_electrodes';

% moviePath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_movie';
% patternsPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_patterns';
% electrodesPath = '/Users/lhruby/MATLAB_code/stimuli/stimulus_files_new/moving_bar_electrodes';
 
%moviePath = '2sec_cluster_movie';
%patternsPath = '2sec_cluster_patterns';
%electrodesPath = '2sec_cluster_electrodes';

%moviePath = '1el_local_movie';
%patternsPath = '1el_local_patterns';
%electrodesPath = '1el_local_electrodes';

%moviePath = '/Users/lhruby/Documents/stim_files_backup/2011-01-11-0/inputs/2sec_cluster_movie';
%patternsPath = '/Users/lhruby/Documents/stim_files_backup/2011-01-11-0/inputs/2sec_cluster_patterns';
%electrodesPath = '/Users/lhruby/Documents/stim_files_backup/2011-01-11-0/inputs/2sec_cluster_electrodes';

moviePath = 'spatiotemp_probe_movie';
patternsPath = 'spatiotemp_probe_patterns';
electrodesPath = 'spatiotemp_probe_electrodes';

%moviePath = '/Users/lhruby/Documents/stim_files_backup/2010-05-14/2sec_cluster_movie';
%patternsPath = '/Users/lhruby/Documents/stim_files_backup/2010-05-14/2sec_cluster_patterns';
%electrodesPath = '/Users/lhruby/Documents/stim_files_backup/2010-05-14/2sec_cluster_electrodes';

%moviePath = '/Users/lhruby/Documents/stim_files_backup/2011-07-14-0/inputs/moving_bar_movie';
%patternsPath = '/Users/lhruby/Documents/stim_files_backup/2011-07-14-0/inputs/moving_bar_patterns';
%electrodesPath = '/Users/lhruby/Documents/stim_files_backup/2011-07-14-0/inputs/moving_bar_electrodes';

 
updatePeriod = 20; % (in samples)


%% decoding electrodes file and patterns file

electrodesID = fopen(electrodesPath, 'r', 'ieee-le.l64'); %read permission, little-endian
electrodes = fread(electrodesID, inf, 'int32');
fclose(electrodesID);

nElec = length(electrodes);

patternsID = fopen(patternsPath, 'r', 'ieee-le.l64'); %read permission, little-endian
patterns = fread(patternsID, inf, 'double');
fclose(patternsID);

patterns = reshape(patterns, nElec, []);

maxAmp = max(max(abs(patterns)));

patternsNorm = patterns/maxAmp;

keyboard

%% decoding movie chunks file

movieID = fopen(moviePath, 'r', 'ieee-le.l64'); %read permission, little-endian
movieChunksFile = fread(movieID, inf, 'int32');
fclose(movieID);

nChunks = movieChunksFile(1);
movieChunks = cell(nChunks, 1);
chunkSamples = zeros(nChunks, 1);

startIndex = 2;
for i = 1:nChunks
    chunkSize = movieChunksFile(startIndex)-6; %first six samples are different
    chunkSamples(i) = movieChunksFile(startIndex+6);
    startIndex = startIndex+7; %start of (time, pattern, gain) values
    movieChunks{i} = zeros(chunkSize/3, 2); %first column = times, second column = patterns
    movieChunks{i}(:,1) = movieChunksFile(startIndex:3:startIndex+chunkSize-3);
    movieChunks{i}(:,2) = movieChunksFile(startIndex+1:3:startIndex+chunkSize-2);
    %for now, gain values are ignored (not used)
    startIndex = startIndex + chunkSize; %start of next movie chunk
end

keyboard

%% plot as rasters

if 0
    nPatterns = length(patterns);
    chosenChunk = 4;
    
    figure
    for i = 1:nPatterns
        axes('position', [0.1 (nPatterns - i + 1)/(nPatterns+1) 0.8 1/(nPatterns+2)])
        hold on
        thesePatternInd = find(movieChunks{chosenChunk}(:,2)==i);
        for j = 1:length(thesePatternInd)
            plot([movieChunks{chosenChunk}(thesePatternInd(j),1) movieChunks{chosenChunk}(thesePatternInd(j),1)]/20000, [0 1], 'k-', 'LineWidth', 1)
        end
        hold off

        set(gca, 'yLim', [0 1], 'xlim', [0 1])
        if i == nPatterns
            xlabel('time (s)')
        else
            set(gca, 'xtick', [])
        end
    end

end

%%

movieChunksSparse = cell(nChunks, 1);

for i = 1:nChunks
    movieChunksSparse{i} = cell(ceil(chunkSamples(i)/updatePeriod)+1, 1);
    for j = 1:length(movieChunks{i})
        iFrame = floor(movieChunks{i}(j,1)/updatePeriod)+1;
        movieChunksSparse{i}{iFrame} = [movieChunksSparse{i}{iFrame} movieChunks{i}(j,2)];
    end
end


%% generate movie
[xCoords yCoords] = getElectrodeCoords61();
cathColor = [50 70 247]/255;
anColor = [.8 .05 0.05];

stimedSoFar = zeros(64,1);
movieFrames = cell(nChunks, 1);
for i = 1:nChunks
    figure('position', [50 50 600 600])
    movieAxes = axes();
    for j = 1:length(movieChunksSparse{i})
        axes(movieAxes)
        cla
        hold on
        plot(xCoords, yCoords, 'ko', 'MarkerSize', 5)
        plot(xCoords(mod(stimedSoFar,4)==1), yCoords(mod(stimedSoFar,4)==1), 'ro', 'MarkerSize', 5)
        plot(xCoords(mod(stimedSoFar,4)==2), yCoords(mod(stimedSoFar,4)==2), 'mo', 'MarkerSize', 5)
        plot(xCoords(mod(stimedSoFar,4)==3), yCoords(mod(stimedSoFar,4)==3), 'bo', 'MarkerSize', 5)
        plot(xCoords(~~((mod(stimedSoFar,4)==0) .* (stimedSoFar~=0))), yCoords(~~((mod(stimedSoFar,4)==0) .* (stimedSoFar~=0))), 'go', 'MarkerSize', 5)
        fastFrame = 1;
        for k = 1:length(movieChunksSparse{i}{j})
            for m = 1:nElec
                if patternsNorm(m, movieChunksSparse{i}{j}(k)) > 0 %cathodal
                    plot(xCoords(electrodes(m)), yCoords(electrodes(m)), 'o',...
                        'MarkerSize', 50*abs(patternsNorm(m, movieChunksSparse{i}{j}(k))),...
                        'MarkerEdgeColor', cathColor, 'MarkerFaceColor', cathColor)
                    fastFrame = 0;
                    stimedSoFar(electrodes(m)) = stimedSoFar(electrodes(m)) + 1;
                elseif patternsNorm(m, movieChunksSparse{i}{j}(k)) < 0 %anodal
                    plot(xCoords(electrodes(m)), yCoords(electrodes(m)), 'o',...
                        'MarkerSize', 50*abs(patternsNorm(m, movieChunksSparse{i}{j}(k))),...
                        'MarkerEdgeColor', anColor, 'MarkerFaceColor', anColor)
                    fastFrame = 0;
                end
            end
        end
        title(['frame ' num2str(j-1) 10 '(' num2str((j-1)*updatePeriod/20) ' ms)'])
        hold off
        set(movieAxes, 'xlim', [-10 10], 'ylim', [-10 10])
        axis equal
        if fastFrame
            pause(0.05)
        else
            pause(0.1)
        end
        %movieFrames{i} = getframe;
        %M(j) = getframe;
    end
    %movie(M)
end

