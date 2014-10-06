function spectrum=NS_PlotMultipleSpectrumsOnFigure(Traces,FigureProperties,NS_GlobalConstants);
ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

Colors=FigureProperties.Colors;

[f,spectrum]=NS_CalculateMultipleFFTs(Traces,Fs);

s=size(spectrum);
for i=1:s(1)
    signal=abs(spectrum(i,:));
    %signal=Traces(i,:);
    WaveformColorIndex=i-floor(i/length(Colors))*length(Colors)+1;    
    semilogy(f,signal,Colors(WaveformColorIndex));  
    hold on;
end