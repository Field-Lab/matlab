function [CorrectedTraces,EI,Noise]=NS512_CorrectSignaturesForTiming(TracesWithSpikes,SpikesIndexes);
%UNTITLED Summary of this function goes here
%   Funkcja zwraca macierz przebieg�w o zuniformowanych w dziedzinie czasu.
%   Dokonuje si? korekcji czasowej macierzy channelTraces2D na podstawie 
%   wskaznikow podanych w macierzy UniSpikesIndicator  
% ChannelTraces2D - wejsciowa macierz 2-D amplituda vs. pr�bka
% UniformSpikes - wyjsciowa macierz Spikeow uporzadkowanych
% Noise - wyjsciowa macierz przebiegow nie zakwalifikowanych jako spike-i
% Mean - wektor zawierajacy przebieg usredniony 

[NumberOfTraces NumberOfChannels NumberOfSamples]=size(TracesWithSpikes);
L=25;

if NumberOfSamples-max(SpikesIndexes)<=L
    CorrectedTraces=zeros(NumberOfTraces,NumberOfChannels,L);
    EI=zeros(NumberOfChannels,L);
    Noise=std(CorrectedTraces,0,3);
    warning('The spike delay is too large');
    return;
end

for i=1:NumberOfTraces %Przemiatanie po wszystkich przebiegach
    ind=SpikesIndexes(i); %which sample is first for given pulse        
    if(NumberOfSamples-ind>L)
        %CorrectedTraces(i,:,:)=TracesWithSpikes(i,:,NumberOfSamples-L+1:NumberOfSamples);                
        CorrectedTraces(i,:,:)=TracesWithSpikes(i,:,ind:ind+L-1);
    else
        %error(['The spike delay is too large']);
    end
end

SCT=size(CorrectedTraces);
EI=reshape(mean(CorrectedTraces),SCT(2),SCT(3));
Noise=std(CorrectedTraces,0,3);