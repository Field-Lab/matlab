function NewElectrodes=NS512_519InverseMess(Electrodes);
a=importdata('C:\pawel\nauka\granty\NCN_Harmonia2013\519mess.txt');
for i=1:length(Electrodes)
    NewElectrodes(i)=find(a==Electrodes(i));
end