cd D:\2008-08-26-0\;
SPfilename='pattern006';
PatternNumber=19;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);

Parasol50Movies=[34:3:58];
Parasol50Eff=[ 5 10 14 18 34 49 66 87 95];

Parasol50Amp=[0.5813 0.6319 0.6824 0.7583 0.8341 0.9099 1.0110 1.1043 1.2047];
for i=1:0
    Movie=Parasol50Movies(i);
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,Movie,NS_GlobalConstants);
    [Name,Channels,AmpOut]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants);
    Parasol50Amp(i)=AmpOut;
end
figure(10);
plot(Parasol50Amp,Parasol50Eff,'bd-');

Parasol100Eff=[6 11 14 20 21 34 47 58 68 75 92 99];
Parasol100Movies=[26:3:59];
Parasol100Amp=[0.2136 0.2388 0.2576 0.2780 0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066];
for i=1:0 %length(Parasol100Movies)
    Movie=Parasol100Movies(i);
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,Movie,NS_GlobalConstants);
    [Name,Channels,AmpOut]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants);
    Parasol100Amp(i)=AmpOut;
end
figure(10);
plot(Parasol50Amp/2,Parasol50Eff,'bd-',Parasol100Amp,Parasol100Eff,'rd-');

Midget50Movies=[37:3:61];
Midget50Eff=[4 13 25 31 34 44 74 85 89];

Midget50Amp=[0.6319 0.6824 0.7583 0.8341 0.9099 1.0110 1.1043 1.2047 1.3051];
for i=1:0 %length(Midget50Movies)
    Movie=Midget50Movies(i);
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,Movie,NS_GlobalConstants);
    [Name,Channels,AmpOut]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants);
    Midget50Amp(i)=AmpOut;
end 
plot(Parasol50Amp/2,Parasol50Eff,'bd-',Parasol100Amp,Parasol100Eff,'rd-',Midget50Amp/2,Midget50Eff,'gd-');

Midget100Eff=[11 18 23 34 38 47 59 71 90 95];
Midget100Movies=[32:3:59];
Midget100Amp=[0.2576 0.2780 0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066];
for i=1:0 %length(Parasol100Movies)
    Movie=Parasol100Movies(i);
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,Movie,NS_GlobalConstants);
    [Name,Channels,AmpOut]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants);
    Parasol100Amp(i)=AmpOut;
end
h=plot(Parasol50Amp,Parasol50Eff,'bd-',Parasol100Amp,Parasol100Eff,'bd--',Midget50Amp,Midget50Eff,'rd-',Midget100Amp,Midget100Eff,'rd--');
set(h,'LineWidth',2);
grid on;
FS=16;
h=gca;
set(h,'FontSize',FS);
h=xlabel('Amplitude [\muA]');
set(h,'FontSize',FS);
h=ylabel('Efficacy [%]');
set(h,'FontSize',FS);

figure(11);
h=plot(Parasol50Amp*0.05,Parasol50Eff,'bd-',Parasol100Amp*0.1,Parasol100Eff,'bd--',Midget50Amp*0.05,Midget50Eff,'rd-',Midget100Amp*0.1,Midget100Eff,'rd--');
set(h,'LineWidth',2);
grid on;
FS=16;
h=gca;
set(h,'FontSize',FS);
h=xlabel('Amplitude [nC]');
set(h,'FontSize',FS);
h=ylabel('Efficacy [%]');
set(h,'FontSize',FS);
