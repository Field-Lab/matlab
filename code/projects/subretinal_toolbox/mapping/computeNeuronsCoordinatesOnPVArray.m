function [neuronPos neuronIDList] = computeNeuronsCoordinatesOnPVArray(visionDataFolder, MEAtoArray, neuronIDList)
% [neuronPos neuronIDList] = computeNeuronsCoordinatesOnPVArray(visionDataFolder,...
%                               MEAtoArray, neuronIDList);
%
% This function computes the coordinates of the neurons as estimated from
% the EI in the array coordinates system. 
%
% Parameters:
%   - visionDataFolder: path to the folder in which the .params file for
%   the dataset can be found, as well as the .neurons file if the
%   neuronIDList is not specified.
%   - MEAtoArray: coordinates transformation matrix as returned by the
%   mapCoordinates function.
%   - neuronIDList: optional, a list of the neurons for which the new
%   coordinates have to be computed. If it is not specified, all the
%   neurons with valid EI coordinates in the parameter file will be used.
%
% Returns:
%   - neuronPos: a nNeurons*2 matrix containing the coordinates of all the
%   neurons in the pixel array system.
%   - neuronIDList: a nNeuron vector of neuron IDs, so that the neuron with
%   coordinates neuronPos(kk,:) has the ID neuronIDList(kk).
%

% Version: 1.0 - 05/31/2013
% Author: Georges Goetz - Stanford University

%% Reading the input arguments

% Making sure folders ends by '\' or '/'
if visionDataFolder(end:end)~=filesep
    visionDataFolder = [visionDataFolder filesep];
end

%% Linking to the params and neurons files

if isunix
    javaaddpath '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    javaaddpath 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end
if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

% Finding the neuron file, the param file and the ei file
contentsDataFolder = dir(visionDataFolder);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isParamFile = strfind(contentsDataFolder(kk).name,'.params');
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuron_path =  [visionDataFolder contentsDataFolder(kk).name];
            end
        end
    end
    if isParamFile
        param_path = [visionDataFolder contentsDataFolder(kk).name];
    end
end
% Linking to the files
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
paramFile=edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

%% Getting the old neurons coordinates

if ~exist('neuronIDList','var')||isempty(neuronIDList)
    neuronIDList = neuronFile.getIDList();
end
nNeurons = length(neuronIDList);

neuronsMEACoordinates = zeros(nNeurons,2);
for kk=1:nNeurons
    neuronsMEACoordinates(kk,1)=paramFile.getDoubleCell(neuronIDList(kk), 'EIx0');
    neuronsMEACoordinates(kk,2)=paramFile.getDoubleCell(neuronIDList(kk), 'EIy0');
end

%% Transforming them into the pixel array coordinates

neuronPos = zeros(size(neuronsMEACoordinates));
for kk=1:nNeurons
    neuronPos(kk,:) = (MEAtoArray.R*neuronsMEACoordinates(kk,:).' + MEAtoArray.t).';
end

end % computeNeuronsCoordinatesOnArray