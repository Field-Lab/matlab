function EI=NS512Read_EI_File(FilePath,NeuronID,Channels);

fid2=fopen(FilePath,'r');
a=fread(fid2,'integer*2');
fclose(fid2);

NeuronIDs=a(4:3+a(3));
WhichNeuron=find(NeuronIDs==NeuronID)
Index=3+a(3)+(WhichNeuron-1)*a(1)*a(2)+1;
EIfull=reshape(a(Index:Index+a(1)*a(2)-1),a(1),a(2));
EI=EIfull(Channels,:);