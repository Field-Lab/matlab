%% ideas

%% preprocessing
cd('/mnt/muench_data/user/alexandra/scripts/spike_sorting')
clc
close all
clear all
clear functions
clear java
clear classes
pack
pause(1)

path(path,'/mnt/muench_data/user/alexandra/scripts/spike_sorting')
path(path,'/mnt/muench_data/user/alexandra/scripts/spike_sorting/InsidePolyFolder')

global spikes spikeTimes hekaNames thr thrs src traceUI getXdata getYdata qualityUI notesUI autoPathUI...
    fileName pathName savePath date ch experimenter

get_file(1)