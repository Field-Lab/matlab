function [UniformSpikes, MeanSpike, Noise] = SpikeUnif_PR(ChannelTraces2D, UniSpikesIndicator, TimeCorelSpikes);
%UNTITLED Summary of this function goes here
%   Funkcja zwraca macierz przebiegów o zuniformowanych w dziedzinie czasu.
%   Dokonuje si? korekcji czasowej macierzy channelTraces2D na podstawie 
%   wskaznikow podanych w macierzy UniSpikesIndicator  
% ChannelTraces2D - wejsciowa macierz 2-D amplituda vs. próbka
% UniformSpikes - wyjsciowa macierz Spikeow uporzadkowanych
% Noise - wyjsciowa macierz przebiegow nie zakwalifikowanych jako spike-i
% Mean - wektor zawierajacy przebieg usredniony 

[NumberOfTraces NumberOfSamples] = size(ChannelTraces2D);

%Odnajdywanie przebiegow
L=34;   %Liczba probek do przepisania na wyjscie

k = 1;  %Pomocnicza zmienna do indeksowania macierzy UniformSpikes
n = 1;  %Pomocnicza zmienna do indeksowania macierzy Noise
%a=find(TimeCorelSpikes == 1);
ChannelTracesUni = zeros(length(find(TimeCorelSpikes==1)), L);
%ChannelTracesUni = zeros(size(find(sign(UniSpikesIndicator)) == 1), 24);
for i=1:NumberOfTraces %Przemiatanie po wszystkich przebiegach
    ind=UniSpikesIndicator(i);
    typ=TimeCorelSpikes(i);
    %UniSpikesIndicator(i)
    if typ==1
        if(NumberOfSamples - ind <L)
            NumberOfSamples - ind;
            ChannelTracesUni(k,:) = ChannelTraces2D(i,NumberOfSamples-L+1:NumberOfSamples);
            k=k+1;
        else 
        ChannelTracesUni(k,:) = ChannelTraces2D(i,ind:ind+L-1);
        k = k+1;
        end
    else
        Noise(n,:) = ChannelTraces2D(i,:);
        n = n+1;
    end
end

 %display('Rozmiar macierzy Uni');
size(ChannelTracesUni);
UniformSpikes = ChannelTracesUni;  % + SpikeThresh;
MeanSpike = mean(UniformSpikes,1);
end























































































































































































