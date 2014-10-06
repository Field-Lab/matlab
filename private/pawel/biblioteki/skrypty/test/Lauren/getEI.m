%%%script to extract maximum voltage amplitudes for each neuron/electrode combination from the ei file

close all;
clear all;

pathToEi = 'E:\analysis\2008-08-26-0\data005\data005.ei';
NeuronID=31;

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ids = eiFile.getIDList()
EI=eiFile.getImage(NeuronID);

s=reshape(EI(1,:,:),65,81);
plot(s')

break;

%pathToNeurons = '/netapp/snle/home/lhruby/Desktop/data001.classification.txt'; %path to a text file containing the ids of the neurons to be analyzed, separated by returns

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ids = eiFile.getIDList();

neurons = load(pathToNeurons);
nsize = size(neurons);
n = squeeze(nsize(:,1)); %number of neurons specified

eiAmps = zeros(n,512); %creates an array of zeros to store the maximum voltage for each neuron/electrode average voltage trace

tic
for i = 1:1:n,    
    for j = 1:1:513,
        tempID = neurons(i); %gets the id number of the ith neuron
        ei = eiFile.getImage(tempID); %gets ei data for ith neuron, storing the information as a 3D array
        trace = squeeze(ei(1,j,:)); %gets average trace of neuron i on electrode j-1 and stores it as a 1D array "trace"
        eiAmps(i,j) = max(trace); %finds the maximum value in the trace and stores it in eiAmps
    end
end
toc

figure(1)
surf(eiAmps)