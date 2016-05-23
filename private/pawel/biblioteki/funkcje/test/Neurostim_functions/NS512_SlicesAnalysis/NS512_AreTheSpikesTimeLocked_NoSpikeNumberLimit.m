function TimeLocking=NS512_AreTheSpikesTimeLocked_NoSpikeNumberLimit(SpikeTimesHistogram,TimeWindow);
% Ta funkcja sprawdza, czy spiki na danej elektrodzie sa czasowo
% skorelowane ze stymulacja. Na wejsciu podawany jest histogram opoznien
% spikoe (Post Stimulus Time Hostogram) okreslony np. dla 600 wartosci
% opoznien. Histogramy sa sumowane w przesuwajacym sie oknie czasowym o
% zadanej dlugosci, i jesli dla pewnego polozenia okna czasowego ilosc
% mieszcacych sie w nim spikow jest co najmniej trzy razy wieksza niz by to
% wynikalo z czystej statystyki (opartej na zalozeniu "przypadkowych"
% momentow pojawienia sie spikow), oraz gdy taka liczba spikow w tym oknie
% jest rowna co najmniej 10, to wtedy zakladamy ze spiki sa skorelowane
% czasowo.
% SpikeTimesHistogram - o wymiarach ilosc_kanalow x ilosc_probek, np. 512 x
% 600
% Time Window - w ja
for i=1:length(SpikeTimesHistogram)-TimeWindow+1
    NumberOfSpikesInWindow(i)=sum(SpikeTimesHistogram(i:i+TimeWindow-1));
end

%sum(SpikeTimesHistogram)
MeanSpikeRate=sum(SpikeTimesHistogram)/length(SpikeTimesHistogram); %average number of spikes per unit time
MeanSpikeNumberPerWindow=MeanSpikeRate*TimeWindow;

if max(NumberOfSpikesInWindow)>5*MeanSpikeNumberPerWindow & max(NumberOfSpikesInWindow)>25
    TimeLocking=1;
else
    TimeLocking=0;
end

%figure(100)
%plot(SpikeTimesHistogram)