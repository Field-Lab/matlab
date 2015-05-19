function selectNeuronsFromRawFile(dataFolder, idList)
% This function creates a neuron file by selecting the neurons whose ID is
% specified in a neurons-raw file. 
% Warning: if there is a .neurons file in the folder, the function will not
% run.
% 
% Parameters: 
%   - dataFolder: folder containing the neurons-raw file
%   - idList: list of neurons to select. Neurons missing from the raw file
%   will be ignored. 
%
% Version:
%

%% Reading the input arguments

% Making sure dataFolder and outputPath end by '\' or '/'
if dataFolder(end:end)~=filesep
    dataFolder = [dataFolder filesep];
end

if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end

%% Linking to the files

if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

% Finding the old neuron file and the prj file
contentsDataFolder = dir(dataFolder);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    if isNeuronFile
        if isNeuronRawFile
            neuron_path =  [dataFolder contentsDataFolder(kk).name];  
        else  
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                err = MException('MATLAB:InputOutputError',...
                    'A .neurons file already exists in the folder');
                throw(err);
            end
        end
    end
end
oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);

% Getting the path to the new neuron file
sep_pos = strfind(dataFolder,filesep);
neuronName = dataFolder((sep_pos(end-1)+1):(sep_pos(end)-1));
newNeuronPath = [dataFolder neuronName '.neurons'];

% Creating the new neuron file
newNeuronFile = edu.ucsc.neurobiology.vision.matlab.ReadNeurons(newNeuronPath, oldNeuronFile);

%% Filling in the new neuron file

idList = intersect(idList,oldNeuronFile.getIDList);
for kk=1:length(idList)
    cNeuronSpikeTimes = oldNeuronFile.getSpikeTimes(idList(kk));
    cNeuronNSpikes = length(cNeuronSpikeTimes);
    cNeuronElectrode = oldNeuronFile.getElectrode(idList(kk));
    
    res = newNeuronFile.addNeuron(cNeuronElectrode, cNeuronSpikeTimes,...
            idList(kk), cNeuronNSpikes);
    if ~res
        err =  MException('VISION:ReadNeurons',...
            'Error writing neuron to file');
        throw(err);
    end
end

oldNeuronFile.close();
newNeuronFile.closeNeurons();

end % selectNeuronsFromRawFile