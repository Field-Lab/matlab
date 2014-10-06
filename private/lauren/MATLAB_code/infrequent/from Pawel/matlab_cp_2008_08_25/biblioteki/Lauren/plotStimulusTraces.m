function [patternTimes,traces] = plotStimulusTraces(thisPatternLong, status, channels, TimeRange, NS_GlobalConstants);

% generates x and y vectors for plotting the stimulus specified by filename, number_of_PD_chunk, and
% patternNumberToPlot
%
% arguments
%   filename: path from current directory to pattern definitions file
%   number_of_PD_chunk: number of chunk within pattern definitions file
%   patternNumberToPlot: number of pattern within pattern definition chunk
%   channels: vector of electrode numbers to be plotted (must include all electrodes used in
%   specified pattern)
%   ArrayID: 1 (for normal 61 arrays)
%   FigureProperties: see NS_PlotClustersOfSignaturesOnArrayLayout
%   NS_GlobalConstants
%
% written by Lauren Hruby, 2008-09-19

Fs=NS_GlobalConstants.SamplingFrequency;

%[patterns_out, patternsIndeces, status] = ReadPatternDataChunk(filename, number_of_PD_chunk, NS_GlobalConstants);

%remove elements of thisPattern in which no stimulation occurs
thisPattern = thisPatternLong;
for i = length(thisPatternLong):-1:1 %looped backwards so that removal of elements works correctly
    if min(thisPatternLong(i).data(1,:) == 0) == 1; %if all DAC values for this electrode = 0
        thisPattern(i) = [];
    end
end

elecInPattern = length(thisPattern);

%get the current ranges for each electrode included in the pattern
currentRanges = NS_GlobalConstants.CurrentRanges;

elecCurrentStep = zeros(1, elecInPattern);
for i = 1:elecInPattern
    elecCurrentStep(i) = currentRanges(status.ChannelsStatus(thisPattern(i).channel).range + 1)/127;
end

%get channel numbers from thisPattern and finds where they are within argument "channels"
channelsWithStim = zeros(1, elecInPattern);
channelsWithStimIndeces = zeros(1, elecInPattern);
for i = 1:elecInPattern
    channelsWithStim(i) = thisPattern(i).channel;
    if max(channels == channelsWithStim(i))<1
        error('Argument "channels" doesn''t include all electrodes use for stimulation.  Aborting plot.')
    end
    channelsWithStimIndeces(i) = find(channels == channelsWithStim(i));
end

%create vector representing current pulse on each channel to be plotted

patternLength = size(thisPattern(1).data,2);

patternVectors = zeros(length(channels), patternLength*2+1);
patternTimes = zeros(1, patternLength*2+1);

for i = 1:elecInPattern
    patternVectors(channelsWithStimIndeces(i),1) = thisPattern(i).data(1,1)*thisPattern(i).data(3,1)*elecCurrentStep(i);
    patternTimes(1) = 0;
    for j = 1:patternLength-1
        patternVectors(channelsWithStimIndeces(i),j*2) = thisPattern(i).data(1,j)*thisPattern(i).data(3,j)*elecCurrentStep(i);
        patternTimes(j*2) = j;
        patternVectors(channelsWithStimIndeces(i),j*2+1) = thisPattern(i).data(1,j+1)*thisPattern(i).data(3,j+1)*elecCurrentStep(i);
        patternTimes(j*2+1) = j;
    end
    patternVectors(channelsWithStimIndeces(i),patternLength*2) = thisPattern(i).data(1,patternLength)*thisPattern(i).data(3,patternLength)*elecCurrentStep(i);
    patternTimes(patternLength*2) = patternLength;
    %so that vector is extended for entire time window to be plotted
    patternTimes(patternLength*2+1) = TimeRange(2);
end
patternTimes = patternTimes./Fs*1000;

traces = zeros(1,size(patternVectors,1),size(patternVectors,2));
traces(1,:,:) = patternVectors;

%figHandle = NS_PlotClustersOfSignaturesOnArrayLayout(traces,channels,0,ArrayID,FigureProperties,NS_GlobalConstants,patternTimes);