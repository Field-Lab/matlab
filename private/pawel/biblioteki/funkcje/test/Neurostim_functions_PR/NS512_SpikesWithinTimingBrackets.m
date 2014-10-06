function [Types,Indexes]=NS512_SpikesWithinTimingBrackets(SpikesTimings,stdDevPrev,stdDevFoll);

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   stdDevprev - ile odchylen standardowych poprzedzajacych sredni czas
%   detekcji przekroczenia progu
%   stdDevFoll - ile odchylen standardowych nastepujacych po srednim czasie
%   detekcji przekroczenia progu
%   SpikesTimings - wektorzawierajacy wskazniki przekroczenia progu, 0
%   jesli nie nastapilo przekroczenie progu
%   SpikesWithinBrackets - wektor o rozmiarze rownym SpikesTimings,
%   zawierajacy 1 jesli indeks miesci sie w zdefiniowanym zakresie i 0
%   jesli jest poza zakresem stdDevPrev-stdDevFoll

%Znajdywanie niezerowych wskaznikow dla obliczenia sredniej i stdDev
NonZeroUnifiedSpikes = SpikesTimings(SpikesTimings~=0);

if (length(find(SpikesTimings==0))>0)
    warning('Some spikes have timing equal to zero');
end

MeanDelay = mean(NonZeroUnifiedSpikes);
StanDev = std(NonZeroUnifiedSpikes);
Types = SpikesTimings>=(MeanDelay-stdDevPrev*StanDev) & SpikesTimings<=(MeanDelay+stdDevFoll*StanDev);
Indexes=find(Types==1);