function NewElectrodes=NS512_519ApplyMess(Electrodes);
a=importdata('C:\home\Pawel\nauka\Stanford\519mess.txt');
%for i=1:length(Electrodes)
%    NewElectrodes(i)=find(a==Electrodes(i));
%end
NewElectrodes=a(Electrodes);