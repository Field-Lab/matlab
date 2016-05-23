clear;
%
%FilesPath='C:\pawel\nauka\Neurostim-3\sim_files'; % laptop
FilesPath='C:\home\Pawel\nauka\Neurostim2\symulacje1'; % pokoj 023

NumberOfChannels=3; %numbered from 1

SamplingPeriod=25e-6;

HoldLength=10;
FrameHeaderLength=10;
TrigerDuration=4;

FullDataStream=[];
FullClkStream=[];
FullTriggerConnectStream=[];
FullTriggerRecordStream=[];

%Cadence_sim_config1;
%Cadence_sim_realtime1;

Cadence_sim_config2;
Cadence_sim_realtime1;

[clk,data]=Neurostim2_Bitstream2VoltageSignals3(FullDataStream,FullClkStream,FullTriggerConnectStream,FullTriggerRecordStream,25e-9,FilesPath,12,50);