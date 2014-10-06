
%% 
% Get pattern number
patt = 358;

% Get movie number
mov = 216;

% cell
cell_no=7507;

%% 
% open raw file
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(sprintf('/Volumes/Analysis/2012-09-24-3/data008/p%d/p%d_m%d',patt,patt,mov));

%gets data from sample 0 to 20000 (first second)
data = rawFile.getData(0, 20000);

%plot TTLs (electrode 0) ?Recall that matlab starts counting from 1, while
%java starts from 0
plot(1:1:20000, data(:,1));

%plot electrode 1
plot(1:1:20000, data(:,2));

%%
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile('/Volumes/Analysis/2012-09-24-3/data007/data007.ei');
ids = eiFile.getIDList();

ei = eiFile.getImage(cell_no,0)  % (id, errorType (Standard Deviation of the Mean: 0, Variance of the Mean: 1))
% for standard version of matlab, no second argument (standard deviation vs. variance not an option)

% returns ei(average:1 error:2, electrode + 1, time index)
% 3 dimensional array with first dimension referring to whether value is the average voltage (1) or the error (2), second dimension referring to electrode number, and % third dimension referring to point in time

maxElectrode = eiFile.getMaxElectrode(ei);

%% 
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(sprintf('/Volumes/Analysis/2012-09-24-3/data008/p%d',patt),[],0,patt,mov,20000);

iTrial=13;
data =reshape( DataTraces(iTrial,maxElectrode,:),[1,100]);
figure;plot(data);


