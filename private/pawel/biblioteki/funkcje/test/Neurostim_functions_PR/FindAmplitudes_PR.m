function [ Amplitudes_out ] = FindAmplitudes_PR( DataPath,pattern,Movies,Channels,NS_GlobalConstants )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:2:max(Movies)

[~,Amplitudes]=NS_StimulatedChannels(DataPath,pattern,i,[1:512],NS_GlobalConstants);
Amplitudes_out((i+1)/2) = Amplitudes(1);
end

