function [ indexes ] = FindIndexes_PR( Channels, ChannelsPlot )
%UNTITLED Summary of this function goes here
%   Funkcja zwraca indeksy do macierzy zawierajacej dane ze "sprawnymi"
%   elektrodami. Jesli dana elektroda zostala odrzucona w jej miejsce
%   zwraca 0
% Channels - wejsciowa macierz aktywnych kanalow
% ChannelsPlot - wejsciowa macierz kanalow do rysowania
% indexes - wyjsciowa macierz indeksow, w ktorej pozycji macierzy Channels
% znajduja sie wartosci z macierzy ChannelsPlot

ChannelsPlot_size = length(ChannelsPlot);
indexes = zeros(ChannelsPlot_size,1);
for i=1:ChannelsPlot_size
    h = find(Channels == ChannelsPlot(i));
    if isempty(h)
    indexes(i) = 0;
    else
    indexes(i) = h;
    end
end

end

