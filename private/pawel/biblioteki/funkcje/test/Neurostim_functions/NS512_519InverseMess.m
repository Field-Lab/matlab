function NewElectrodes=NS512_519InverseMess(Electrodes);
a=importdata('/Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_files/519mess.txt');
for i=1:length(Electrodes)
    NewElectrodes(i)=find(a==Electrodes(i));
end