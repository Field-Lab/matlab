% this file contains self contained examples about how to call vision in order
% to accomplish a variety of standard tasks. in order to try one, simply
% update the paths in a given example, and then right click on that example
% and go to "Evaluate Current Cell"
%
% if you are interested in improving spike sorting, most steps of the spike
% sorting framework can be directly accessed using the methods here.
%
% before reimplementing any of the tasks below, note that there might
% already be a function that accomplishes the same thing in the standard
% matlab framework (/snle/home/snl-e/matlab-standard/code/lab). when an
% example references another function (e.g. see compute_sta_) it is
% referring to an m file in the path above.
%
% should give you 1.5 or better, otherwise upgrade Matlab
% version -java
%
% if you type javaclasspath the following path should be there:
% /snle/lab/Development/RRS/vision-package/vision.app/Contents/Resources/Java/Vision.jar
%
% tamachado (2009-06-19)
% jgauthier (2008-12-10)
% phli 2012-08-27: Moved non-Vision-specific general patterns to GENERAL_JAVA_CALLS

%%
%-------------------------------------------------------------------------
% open raw data file (read a single .bin file)
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile('/Data/Machado/2009-04-13-0/data000/data000.bin');

% get sampling frequency (in Hz)
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% gets data from sample 0 to 20000 (first second) across ALL electrodes
data = rawFile.getData(0, samplingRate);

% plot TTLs (electrode 0),  Recall that matlab starts counting from 1, while java starts from 0
figure; plot(1:1:samplingRate, data(:,1));

% plot electrode 1
figure; plot(1:1:samplingRate, data(:,2));

%%
%-------------------------------------------------------------------------
% open raw data file (alternate method that reads across data files)
rawFile = edu.ucsc.neurobiology.vision.io.RawDataWrapper('/Data/Machado/2009-04-13-0/data000/data000.bin');

% get sampling frequency (in Hz)
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

%set up arguments
sample = 1000; %sample to start reading from
nSamples = 1*samplingRate; % get one second, starting from sample
nElectrodes = 65; % ALL electrodes in the dataset

%get the raw data
data = rawFile.getData(sample, nSamples, nElectrodes);

% plot TTLs (electrode 0),  Recall that matlab starts counting from 1, while java starts from 0
figure; plot(1:1:samplingRate, data(:,1));

% plot electrode 1
figure; plot(1:1:samplingRate, data(:,2));

%%
%-------------------------------------------------------------------------
% make a new model file and add some clusters to it (see imm-cluster in the spike sorting projects folder)

% create a new model file. note that you need to provide another model file 
% that contains the eigenvectors you want stored in the new model file
oldModelFile = '/Data/Machado/2009-04-13-0/data000/test/test.model';
newModelFile = '/Data/Machado/2009-04-13-0/data000/test/new.model';
newModelObject = edu.ucsc.neurobiology.vision.matlab.ReadModel(newModelFile, oldModelFile);

% open the old model file to get some clusters out of it (this is an
% example, you could generate clusters in any way you want).
electrode = 1;
model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile('/Data/Machado/2009-04-13-0/data000/test/test.model');
% get the clusters on electrode 1
clusters = model.getNeuronExtraction(electrode);
nClusters = length(clusters.probability);
newModelObject.addElectrode(clusters.probability, clusters.means, clusters.covariances, electrode, nClusters)
newModelObject.closeModel;

electrode = 1;
model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile('/Data/Machado/2009-04-13-0/data000/test/new.model');
% get the clusters on electrode 1
clusters = model.getNeuronExtraction(electrode);
nClusters = length(clusters.probability);

% plot the output verifying that everything worked
projectionsPath = '/Data/Machado/2009-04-13-0/data000/test/test.prj';
figure; plot_projections(projectionsPath, 'clusters', newModelFile, 'electrode', electrode);

%%
%-------------------------------------------------------------------------
% read in a model file plot the clusters (see imm-cluster in the spike sorting projects folder)

% see examples above and below this one (or look at plot_projections)

%%
%-------------------------------------------------------------------------
% read in a projections file, a model file, and plot clusters (see imm-cluster in the spike sorting projects folder)

% open the projections file
prjObject = edu.ucsc.neurobiology.vision.matlab.ReadProjections('/Data/Machado/2009-04-13-0/data000/test/test.prj');

% go to electrode 1
electrode = 1;
prjObject.readElectrode(electrode);

% get the projections from electrode 1
prj = prjObject.getProjections();
% get the number of spikes on that electrode
nSpikes = prjObject.getSpikeCount();
prjs = prj(:,1:nSpikes);
% open the model file
model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile('/Data/Machado/2009-04-13-0/data000/test/test.model');
% get the clusters on electrode 1
clusters = model.getNeuronExtraction(electrode);
% plot the projections
figure; plot_projections(prjs, 'clusters', clusters);

% plot_projections will also work by passing in paths
% to a projection file and a model file (see plot_projections)

%%
%-------------------------------------------------------------------------
% open a neurons file (to get spike times for individual neurons)
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('/Data/Machado/2009-04-13-0/data000/test/test.neurons');

% get list of neuron numbers
idList = neuronFile.getIDList();

% get spike times for first neuron
spikeTimes = neuronFile.getSpikeTimes(idList(1));

% histogram of spike times
hist(spikeTimes, 0:10000:double(max(spikeTimes)));

%%
%-------------------------------------------------------------------------
% write out a neurons file (to change which spikes belong to which neurons)
% see (duplicate-removal or imm-cluster in the spike sorting folder)

% make the new neurons file. steal the vision header from a projections
% file or another neurons file...
oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('/Data/Machado/2009-04-13-0/data000/test/test.neurons');
newNeuronFile = edu.ucsc.neurobiology.vision.matlab.ReadNeurons('/Data/Machado/2009-04-13-0/data000/test/new.neurons',oldNeuronFile);

% where was the neuron found
electrode = 1;
% it was the nth cluster on electrode
clusterIndex = 0; 
% specify its spike times
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('/Data/Machado/2009-04-13-0/data000/test/test.neurons');
idList = neuronFile.getIDList();
spikeTimes = neuronFile.getSpikeTimes(idList(1));
% how many spikes does it have
nSpikes = length(spikeTimes);

% add the neuron
newNeuronFile.addNeuron(electrode, clusterIndex, spikeTimes, nSpikes);

% load and plot the spike times verifying that everything worked
newNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('/Data/Machado/2009-04-13-0/data000/test/new.neurons');
idList = newNeuronFile.getIDList();
spikeTimes = newNeuronFile.getSpikeTimes(idList(1));
hist(spikeTimes, 0:10000:double(max(spikeTimes)));

%%
%-------------------------------------------------------------------------
% open a spikes file (to get all spike times)
spikeFile = edu.ucsc.neurobiology.vision.io.SpikeFile('/Data/Machado/2009-04-13-0/data000/test/test.spikes');

% get spikeTimes for electrode 1
spikeTimes = spikeFile.getSpikeTimes(1);

% histogram of spike times
hist(spikeTimes, 0:10000:double(max(spikeTimes)));

%%
%-------------------------------------------------------------------------
% open a params file (to get gaussian fit information)
paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile('/Data/Machado/2009-04-13-0/data000/test/test.params');

% gets list of neuron ids
ids = paramsFile.getIDList();

% to get the class name for a particular neuron.  You can loop through all the neurons to get all the class names.
% it would be nicer to get an array of strings directly.  I do not know if matlab could handle this.
name = paramsFile.getStringCell(ids(1), 'classID');

% gets neurons from a particular class
paramsFile.getNeuronsInClass('All/On/Parasol')

% to see what parameters are possible, look in vision.analysis package
% for classes that implement ParametersCalculator.  The parameters will be in
% the getParameterTypes() function.

% get gaussian fit parameters for a given neuron.
% other Double parameters are gotten the same way.
paramsFile.getDoubleCell(ids(1), 'x0')
paramsFile.getDoubleCell(ids(1), 'y0')
paramsFile.getDoubleCell(ids(1), 'SigmaX')
paramsFile.getDoubleCell(ids(1), 'SigmaY')
paramsFile.getDoubleCell(ids(1), 'Theta')

% get red timecourse.  Other DoubleArray Parameters are gotten the same way.
rt = paramsFile.getArrayCell(ids(1), 'RedTimeCourse');
figure; plot(rt,'r');

%%
%-------------------------------------------------------------------------
% load EI File (see load_ei)
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile('/Data/Machado/2009-04-13-0/data000/test/test.ei');
ids = eiFile.getIDList();

% ei(id, errorType (Standard Deviation of the Mean: 0, Variance of the Mean: 1))
% returns ei(average:1 error:2, electrode + 1, time index)
ei = eiFile.getImage(ids(1),0);  

maxElectrode = eiFile.getMaxElectrode(ei);

%%
%-------------------------------------------------------------------------
% compute EIs (see compute_ei_)
datarun = load_data('/Data/Machado/2009-04-13-0/data000/test/test');
datarun = load_neurons(datarun);

% set parameters
dPath = '/Data/Machado/2009-04-13-0/data000/data000.bin';
nlPoints = 20; % samples before spike to include in EI
nrPoints = 60; % samples after spike to include in EI

% get sampling rate
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% compute eis ftw
ei = edu.ucsc.neurobiology.vision.matlab.Matlab.computeElectrophysiologicalImage(dPath, datarun.spikes{1}*samplingRate, nlPoints, nrPoints, length(datarun.spikes{1}));

%%
%-------------------------------------------------------------------------
% to get electrode maps
rawArrayID = 1501;
electrode = 12;
radius = 1;

electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(rawArrayID);

% if(arrayID < 500), give 64 map
% if(array >= 500 && arrayID < 1500) give 512 map
% if(array >= 1500 && < 2500) give 519 map

xPoint = electrodeMap.getXPosition(electrode);
yPoint = electrodeMap.getYPosition(electrode);

neighbors = electrodeMap.getAdjacentsTo(electrode, radius);  %radius is in nearest neighbor number (integer) 

%%
%-------------------------------------------------------------------------
% compute STAs and a vision movie (see compute_sta_.m and compute_vision_movie.m)
% get triggers and spike times (using the matlab framework)
datarun = load_data('/Data/Machado/2009-04-13-0/data000/test/test');
datarun = load_neurons(datarun);

% get sampling rate
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% compute vision movie
movie = edu.ucsc.neurobiology.vision.matlab.Matlab.computeMovie('/snle/acquisition/movie-xml/RGB-10-8-0.48-11111.xml', datarun.triggers*samplingRate);

% get movie frame
STAFrame = movie.getFrame(frame - 1);
stimulusFrame = permute(reshape(STAFrame.getBuffer,3,fieldWidth,fieldHeight),[3 2 1]) - .5;

% calculate STA
staDepth = 30;
staOffset = -1; % -1 should mean the spike occurs in the last frame of the STA.  0 would give an extra frame past the spike time at the end.
sta = edu.ucsc.neurobiology.vision.matlab.Matlab.computeSTA(datarun.spikes{1}*samplingRate, datarun.triggers(1)*samplingRate, movie, staDepth, staOffset);

%%
%-------------------------------------------------------------------------
% calculate average frame refresh rate
% get triggers and spike times (using the matlab framework)
datarun = load_data('/Data/Machado/2009-04-13-0/data000/test/test');
datarun = load_neurons(datarun);

% get the average frame refresh rate
edu.ucsc.neurobiology.vision.matlab.Matlab.computeRigFrameTime(datarun.triggers*samplingRate)


%%
%-------------------------------------------------------------------------
% Globals file access
globalsFile = edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile('/Volumes/Rat/Data/Max/2009-04-13-7/data002/data002.globals', 0);
movieParams = globalsFile.getRunTimeMovieParams();
movieParams.pixelsPerStixelX  %etcetera
%	public class RunTimeMovieParams{
%		private static final int TAG = 1;
%		public int pixelsPerStixelX, pixelsPerStixelY;
%		public double width, height; //stixels
%		public double micronsPerStixelX, micronsPerStixelY; 
%		public double xOffset, yOffset; //microns, from bottom left
%		public int interval;
%		public double refreshPeriod; //ms
%		public int nFramesRequired;
%		public int[] droppedFrames;
%       }
%
imageParams = globalsFile.getImageCalibrationParams();
imageParams.micronsPerPixelX %etc
%	public class ImageCalibrationParams {
%		private static final int TAG = 0;
%		public double micronsPerPixelX, micronsPerPixelY;
%		public double centerX, centerY; //microns
%		public boolean flipX, flipY; 
%		public double angle;//radians
%		public int arrayID, arrayPart, arrayNParts;
%       }
