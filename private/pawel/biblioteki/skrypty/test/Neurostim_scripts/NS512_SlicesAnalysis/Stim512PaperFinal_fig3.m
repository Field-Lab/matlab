% Detekcja spikow na potrzeby tych rysunkow jest bardzo uproszczona: nie
% wykrywamy spikow, tylko patrzymy, czy sygnal w danym przebiegu (30 ms)
% przekroczyl prog ujemny 40 oraz prog dodatni 20. to mzoe byc nawet jakis
% dziwny przebieg typu LFP, niewazne. Patrzymy tylko na ilu z 50 przebiegow
% mamy tego typy przekroczenie progu. Prog 40 odpowiada 150 uV.
% Od danych zosta? odj?ty zrtefakt z TTX, a 0.5 ms po impulsie zosta?o
% wyzerowane. Dopiero wtedy by?a detekcja impulsów. KTORE ELEKTRODY
% WYRZUCAMY?? aktualnie stymulujaca oraz najblizsi sasiedzi

clear

MovieFile='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\MovieFiles\2013-12-12-3\movie001';
DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
ArtifactDataPath=DataPath;
%FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2013-12-12-3-PH-2015-01-15';
FigurePath='C:\home\Pawel\nauka\512stim_paper\FiguresProby';

NumberOfAmplitudes=25; %25 for 2013-12-12-3
NumberOfMovieSequences=16; %16 for 2013-12-12-3
NumberOfMovieSequences2=18; %18 for 2013-12-12-3

figure(3)
clf
FontSize=14;
Stim512PaperFinal_fig3A_v4;
Stim512PaperFinal_fig3B_v3;
Stim512PaperFinal_fig3D_v2;

FullName=[FigurePath '\Figure3.tif']; 
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 12]);
set(h,'PaperPosition',[0 0 16 12]); 
print(h, '-dtiff', '-r100', FullName);