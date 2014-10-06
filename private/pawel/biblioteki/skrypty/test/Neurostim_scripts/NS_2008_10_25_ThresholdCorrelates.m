NS_GlobalConstants=NS_GenerateGlobalConstants(61);

CenterChannel=44;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Channels=electrodeMap.getAdjacentsTo(CenterChannel,1)';

%1. Define neuron
pathToEi = 'E:\analysis\2008-08-26-0\data005\data005.ei';
NeuronID=558;
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ids = eiFile.getIDList();
EI0=eiFile.getImage(NeuronID);
EI=reshape(EI0(1,Channels,:),length(Channels),81);
Amplitudes=zeros(1,length(Channels));
for i=1:length(Channels)
    Amplitudes(i)=abs(min(EI(i,:)));
end

%1b) Estimate soma position
[x,y]=NS_EstimateMassCenter(Amplitudes,Channels);
direction=[1 -2 0];

%2. 
DataPath='E:\analysis\2008-08-26-0\data008_proba4';
ClusterFileName='E:\analysis\2008-08-26-0\data008_proba4\clusters008';

Pattern=293;
Movies=[34:5:124];
z=0
Patterns=[293:298 305:310];
for i=1:0 %length(Patterns)
    [I(i),V(i),G(i),Gd(i),D2(i),D2d(i)]=NS_FindFieldPropertiesForThreshold(DataPath,ClusterFileName,Patterns(i),Movies,Channels,x,y,z,direction);
end

In=abs(I/max(abs(I)));
Vn=abs(V/max(abs(V)));
Gn=abs(G/max(abs(G)));
Gdn=abs(Gd/max(abs(Gd)));
D2n=abs(D2/max(abs(D2)));
D2dn=abs(D2d/max(abs(D2d)));

x=[1:12];
plot(x,In,'kd-',x,Vn,'bd-',x,Gn,'rd-') %,x,Gdn,'kd--',x,D2n,'bd--',x,D2dn,'rd--');