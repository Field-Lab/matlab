function plotStimulusMovie(movieChunksFile, Array, electrodes)

movieLength = movieChunksFile(8);

movieCropped = movieChunksFile(9:length(movieChunksFile));
nPulses = length(movieCropped)/3;

[elecCoordsX elecCoordsY] = getElectrodeCoords61();


pulseTimes =     zeros(1, nPulses);
patternNumbers = zeros(1, nPulses);

for i = 1:nPulses
   pulseTimes(i) =     movieCropped(1+(i-1)*3);
   patternNumbers(i) = movieCropped(2+(i-1)*3);
end

arrayBinary = (Array > 0);

%M = zeros(1,movieLength/1000);

for i = 0:500:movieLength
    pulseTimesBinary = (pulseTimes > i).*(pulseTimes <= i+500);
    currentPulses = find(pulseTimesBinary);
    figure(1)
    cla
    hold on
    axis([-8 8 -10 10])
    if ~isempty(currentPulses)
        for j = 1:length(currentPulses)
            currentPatternNumber = patternNumbers(j);
            currentPattern = arrayBinary(currentPatternNumber,:);
            currentElectrodes = electrodes(currentPattern);
            plot(elecCoordsX(currentElectrodes), elecCoordsY(currentElectrodes),'.')
        end
    end
    hold off
    M(i/500 + 1) = getframe(gcf);
end

movie(M,1,20)