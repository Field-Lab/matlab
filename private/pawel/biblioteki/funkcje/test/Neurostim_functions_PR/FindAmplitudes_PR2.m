function [ Amplitudes_out ] = FindAmplitudes_PR2( DataPath,pattern,Movies,Channels,NS_GlobalConstants )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(Movies)

[~,Amplitudes]=NS_StimulatedChannels(DataPath,pattern,Movies(i),[1:512],NS_GlobalConstants);
Amplitudes_out(i) = Amplitudes(1);
end

