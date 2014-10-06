function [ UniSpikesIndicator ] = FindUnifiedSpikes_PR(ChannelTraces2D,SpikeThresh);
%FindUnifiedSpikes Funkcja zwraca wektor zawierajacy numery probek, od ktorych
%zaczyna sie przebieg napiecia, ktory mozna zaklasyfikowac jako spike
%   Detailed explanation goes here
%   UniSpikesIndicator - wektor wskazujacy, ktora wartosc danego przebiegu
%   jako pierwsza przekroczyla prog Spike Thresh, jesli nie przekroczyla to
%   zwraca do odpowiadajacej komorki 0.

ChannelTraces2DSize = size(ChannelTraces2D);

%Ustalenie progu napiecia
%SpikeThresh = -25;

%Odejmowanie progu od przebiegow
ChannelTraces2D = ChannelTraces2D - SpikeThresh;

%znajdywanie, ktory przebieg jest powyzej progu, a ktory ponizej
UniSpikesIndicator = zeros(ChannelTraces2DSize(1,1),1);
for i=1:ChannelTraces2DSize(1,1) %Przemiatanie po wszystkich przebiegach
     for j=1:ChannelTraces2DSize(1,2) %Przemiatanie po warto?ciach
         if sign(ChannelTraces2D(i,j))== -1 || sign(ChannelTraces2D(i,j))==0
             UniSpikesIndicator(i) = j;           
             break
         else
             UniSpikesIndicator(i) = 0;
         end
             
     end
end
    %display('Liczba wykrytych spikow');
    size(find(sign(UniSpikesIndicator)) == 1);
end