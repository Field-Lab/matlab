function [TimeCorelSpikes] = FindTimeCorellatedSpikes_PR(UniSpikesIndicator,stdDevPrev,stdDevFoll);

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   stdDevprev - ile odchylen standardowych poprzedzajacych sredni czas
%   detekcji przekroczenia progu
%   stdDevFoll - ile odchylen standardowych nastepujacych po srednim czasie
%   detekcji przekroczenia progu
%   UniSpikesIndicator - wektorzawierajacy wskazniki przekroczenia progu, 0
%   jesli nie nastapilo przekroczenie progu
%   TimeCorelSpikes - wektor o rozmiarze rownym UniSpikesIndicator,
%   zawierajacy 1 jesli indeks miesci sie w zdefiniowanym zakresie i 0
%   jesli jest poza zakresem stdDevPrev-stdDevFoll

%Znajdywanie niezerowych wskaznikow dla obliczenia sredniej i stdDev
NonZeroUnifiedSpikes = UniSpikesIndicator(UniSpikesIndicator~=0);
MeanDelay = mean(NonZeroUnifiedSpikes)
StanDev = std(NonZeroUnifiedSpikes)
TimeCorelSpikes = UniSpikesIndicator>(MeanDelay-stdDevPrev*StanDev) & UniSpikesIndicator<(MeanDelay+stdDevFoll*StanDev)
figure(109)
plot(TimeCorelSpikes)

for i=1:length(UniSpikesIndicator)
    if UniSpikesIndicator(i) == 0;
        TimeCorelSpikes(i) = 0;
    else if UniSpikesIndicator(i)>=(MeanDelay-stdDevPrev*StanDev) & UniSpikesIndicator(i)<=(MeanDelay+stdDevFoll*StanDev)
        TimeCorelSpikes(i) = 1;
    else TimeCorelSpikes(i) = 2;
        end
    end
end
figure(110)
plot(TimeCorelSpikes)
TimeCorelSpikes = TimeCorelSpikes';
end