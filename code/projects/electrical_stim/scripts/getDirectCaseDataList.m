%This script generates a list of cells above a certain threshold and lists axons/somas around them 

%PARAMETERS--------------------
dataPath = '/Volumes/Analysis/2016-01-05-6/data001/data001';
somaThr = -50; 
axonThr = -20;
load('adj_mat_512.mat'); 
out = {}; %initialize output cell

%Read in data (probably only works on Mac)

datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, 'all');

%{
datarun  = load_data(dataPath); %this isn't doing its job on Windows for some reason, so doing it manually - Sasi
dataPathSplit = strsplit(dataPath, filesep); 
disp(class(dataPathSplit))
dataPathSplit = dataPathSplit(2:(end-1)); %needed for MAC
datarun.names.rrs_neurons_path = [dataPath filesep dataPathSplit{end} '.neurons'];
datarun.names.rrs_ei_path = [dataPath filesep dataPathSplit{end} '.ei'];
datarun  = load_neurons(datarun);
datarun  = load_ei(datarun, 'all');
%}

%get and reshape EI data
numNeurons = size(datarun.ei.eis,1); 
allEIs = zeros(numNeurons, size(datarun.ei.eis{1},1), size(datarun.ei.eis{1},2));
for n = 1:numNeurons
	allEIs(n,:,:) = datarun.ei.eis{n}; 
end
eimat = min(allEIs, [], 3); %get minimum for each neuron at each electrode
[minVals, minElecs] = min(eimat, [], 2); %spike amp and main electrode for each neuron

%find which cells pass threshold & filter eimat
cellIndices = find(minVals < somaThr);
eimat = eimat(cellIndices,:);

%zero out all soma electrodes and two adjacent
for i = 1:length(cellIndices);
	elecsToZero = [];
	mainElec = minElecs(cellIndices(i));
	foo = adj_mat_512{mainElec};
	for j = foo;
		elecsToZero = [elecsToZero adj_mat_512{j}];
	end
	elecsToZero = [elecsToZero foo];
	elecsToZero = [elecsToZero mainElec];
	eimat(i,elecsToZero) = 0;
end

%Loop through all neurons that passed threshold
nvec = 1:length(cellIndices);
cnt = 1;
for n = nvec; 
	mainElec = minElecs(cellIndices(n)); %get principle electrode for this neuron
	elecsToCheck = [adj_mat_512{mainElec} mainElec]; %vector of main electrode and all adjacent 
	foo = eimat(:, elecsToCheck) < axonThr;			
	footmp = max(foo, [], 2); 
	cells = find(footmp);%cell ids that pass criteria
	cells = cells'; 
	if ~isempty(cells);
		elecArr = {};
		cntjr = 1; 
		for c = cells;
			grot = foo(c,:).*elecsToCheck; 
			elecArr{cntjr} = grot(grot~=0); %add vector of relevant electrodes for each involved cell
			cntjr = cntjr + 1;
		end
		out{cnt} = {datarun.cell_ids(cellIndices(n)) datarun.cell_ids(cellIndices(cells)) elecArr}; %output entry in format {cellID [involved cells] {[relevant elecs cell 1] [..cell 2..]..]}}
	cnt = cnt + 1;
	end
end
