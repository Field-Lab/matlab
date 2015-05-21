function newNeuronFile = initVisionNeuronFile(rawDataFilePath, newNeuronFilePath, ttlTimes, visionPath)
%INITVISIONNEURONFILE Initializes a Vision neuron file.
% 
%   INITVISIONNEURONFILE(RAWDATAFILEPATH, NEWNEURONFILEPATH, TTLTIMES, ...
%      VISIONPATH) 
%   abstracts away some of the gritty details associated with
%   creating a neuron file that makes sense to Vision.
%   RAWDATAFILEPATH should be the path to a Vision bin file with a header
%   or a folder of bin files. NEWNEURONFILEPATH is the path where the new
%   neuron file will be created. TTLTIMES is a vector of TTL times as 
%   recorded on channel 0 of the raw data file. 
%   VISIONPATH is optional and can be specified if the Vision jar file 
%   has not yet been added to the java class path.
%
%   The function returns a edu.ucsc.neurobiology.vision.io.NeuronFile object
%   to which neurons can be added using the class method:
%      neuronFile.addNeuron(int electrode, int neuronID, int[] times, 
%          int nSpikes)
%   Note the integer type of the arguments.
%   If you don't know what your neuronID is supposed to be, you can call
%      neuronFile.getNeuronID(electrode, clusterIndex) 
%   and Vision will calculate it for you.
%
%   Don't forget to close() the neuronFile after you're done with it.
%   

% Check we linked to the Vision jar file already
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class') && exist('visionPath', 'var')
    javaaddpath(visionPath);
end

% Get the bin file header 
rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataFilePath);
rawDataHeader = rawDataFile.getHeader();
rawDataFile.close();

% Magic number for the Vision Header. 
% If things break look for the value of 
%   edu.ucsc.neurobiology.vision.util.VisionParams.NEURONS_HEADER_CAPACITY
% in the Vision source code. Not sure how to link to it directly from
% Matlab as there is no public constructor for this class.
NEURONS_HEADER_CAPACITY = 50000;
% Same with the Vision Version value (also found in VisionParams)
VISION_VERSION = 7001002%7003047;%8001010;
% Now the neurons file header version is hidden in:
%   edu.ucsc.neurobiology.vision.io.NeuronFile.INT_VERSION
VERSION = 32;

% Instantiate a new Vision header 
visionHeader = edu.ucsc.neurobiology.vision.io.VisionHeader();

% There's of course a MAGIC number because Vision likes magic.
% As far as I can tell it's supposed to be:
%   edu.ucsc.neurobiology.vision.io.ProjectionFile.MAGIC
% which is currently set to 0xBBBBBB
MAGIC = 12303291%16448250%hex2dec('BBBBBB'); % 

% Then some reasonable defaults. Not sure why they're needed but 
% Vision specifies them and it's not a good idea to not specify what 
% Vision specifies, in general, because then things break in unexpected
% places.
MIN_SPIKES = 100;
MAX_CONTAM = 0.1;

% Fill in the Vision header object
visionHeader.magic = MAGIC;
visionHeader.headerVersion = 1;%1;

% visionHeader.headerCapacity = uint(NEURONS_HEADER_CAPACITY);
visionHeader.version = int32(VERSION);
% visionHeader.meanTimeConstant =0%-1;
% visionHeader.threshold = 0%-1;
% visionHeader.arrayID = rawDataHeader.getArrayID();
visionHeader.nSamples = int32(rawDataHeader.getNumberOfSamples());
visionHeader.samplingFrequency = int32(rawDataHeader.getSamplingFrequency());
visionHeader.visionVersion = int32(VISION_VERSION);
% visionHeader.nDimensions = 0%5;
visionHeader.minNeuronSpikes = double(MIN_SPIKES);
visionHeader.maxContamination = double(MAX_CONTAM);
visionHeader.removeDuplicates = int32(1);%-2;
visionHeader.covarianceType = int32(-2);
% visionHeader.acfT1 = 0 %0.5;
% visionHeader.acfT2 = 0% 1;
% visionHeader.coincidenceTime = 0%1;
% visionHeader.maxCorrelation = 0%1;
% visionHeader.time = long(0)%;
% visionHeader.minCovarianceSpikes =0
% visionHeader.maxCovarianceSpikes =0
% visionHeader.electrodeUsage = 1
% TTL times should be integers. Cast them to int32 just in case
ttlTimes = int32(ttlTimes);

% Now we can instantiate the new neuron file
newNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(newNeuronFilePath, ...
    visionHeader, NEURONS_HEADER_CAPACITY, ttlTimes);

end % initVisionNeuronFile