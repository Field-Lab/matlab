function NewElectrodes = NS512_519Inverse(Electrodes)
a = importdata('519mess.txt');
for ii=1:length(Electrodes)
    NewElectrodes(ii)=find(a==Electrodes(ii));
end