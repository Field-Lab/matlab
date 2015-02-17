prjPath = '/Volumes/Analysis/2014-11-05-3/data000/data000.prj';
prjFile = edu.ucsc.neurobiology.vision.matlab.ReadProjections(prjPath);

prjFile.readElectrode(233);

prjs = prjFile.getProjections();
counts = prjFile.getSpikeCount();
times = prjFile.getSpikeTimes();

figure; scatter(prjs(1,:), prjs(2,:), 1); title('1-2'); %axis([-165 50 -280 50]);
figure; scatter(prjs(1,:), prjs(3,:), 1); title('1-3'); %axis([-165 50 -280 50]);
figure; scatter(prjs(2,:), prjs(3,:), 1); title('2-3'); %axis([-280 40 -280 40]);