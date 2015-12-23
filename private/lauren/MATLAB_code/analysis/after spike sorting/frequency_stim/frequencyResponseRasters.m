function axesH = frequencyResponseRasters(elecRespPath, freq, movieInd, nSequence, analyzedRegion, varargin)

% generates a raster plot of responses to pulse train at particular
% frequency
%
% ****note that 'trials' aren't really meaningful for frequencies that result in
% no gap in time between pulse sequences (e.g. 20 pulses at 20 Hz, repeated
% each second)


p = inputParser;

p.addRequired('elecRespPath', @ischar)
p.addRequired('freq', @isnumeric) %frequency of stimulus in Hz
p.addRequired('movieInd', @isnumeric)
p.addRequired('nSequence', @isnumeric) %number of pulses in a single repetition
p.addRequired('analyzedRegion', @isnumeric) %[start end] of analyzed region relative to stimulus, in samples

p.addParamValue('rasterType', 'line', @(x)any(strcmp(x,{'line','dot'}))) %'line' or 'dot'
p.addParamValue('markerColor', [0 0 0], @isnumeric)
p.addParamValue('markerSize', 2, @isnumeric) %only used for rasterType = dot case

p.parse(elecRespPath, freq, movieInd, nSequence, analyzedRegion, varargin{:})

rasterType = p.Results.rasterType;
markerColor = p.Results.markerColor;
markerSize = p.Results.markerSize;




[latencies successes] = extractFrequencyAnalysis(elecRespPath, movieInd, nSequence);

nReps = size(latencies,2);


period = 20000/freq; %in samples


analyzedRegions = zeros(2,nSequence);
analyzedRegions(1,:) = (0:period:period*(nSequence-1)) + analyzedRegion(1);
analyzedRegions(2,:) = (0:period:period*(nSequence-1)) + analyzedRegion(2);
analyzedRegions = analyzedRegions/20;


%determine response times relative to start of sequence
spikeTimes = zeros(size(latencies));
for ii = 1:nSequence
    spikeTimes(ii,:) = latencies(ii,:) + (ii-1)*period;
end

%re-zero failures
spikeTimes = spikeTimes.*successes;

%convert to ms
spikeTimes = spikeTimes/20;

figure
axesH = axes;
hold on

for jj = 1:nSequence
    p = patch([analyzedRegions(1,jj) analyzedRegions(1,jj) analyzedRegions(2,jj) analyzedRegions(2,jj)],...
        [0.5 nReps+0.5 nReps+0.5 0.5], [0.9 0.9 0.9]);
    set(p, 'edgeColor', 'none')
end

for jj = 1:nReps
    x = spikeTimes(:,jj) ~= 0;
    spikesToPlot = spikeTimes(x,jj);
    if ~isempty(spikesToPlot)
        if strcmpi(rasterType,'dot')
            plot(spikesToPlot, jj, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        else
            for kk = 1:length(spikesToPlot)
                plot([spikesToPlot(kk) spikesToPlot(kk)], [jj-0.4 jj+0.4], 'color', markerColor)
            end
        end
    end
end
hold off
set(gca, 'xLim', [0 period*(nSequence+1)/20], 'yLim', [0 nReps+1], 'ytick', [1:nReps])
xlabel('time (ms)')
ylabel('trial')