function save_neurons(filename, header, triggers, spikes, cell_ids, electrodes)
% saveNeurons  Save out a neurons file
%
% usage: saveNeurons(filename, header, triggers, spikes, cell_ids, electrodes);
%
% shlens 2008-04-17
%greschner change to vision

import('edu.ucsc.neurobiology.vision.matlab.*');
import('edu.ucsc.neurobiology.vision.io.*');

%header size hardcoded in vision
headerSize = 50000;

% update user
fprintf(1,['Saving ' num2str(length(spikes)) ' neurons ... ']);

%HACK
if isempty(triggers)
    triggers=0;
end

% create new neurons file
neuronsFile = NeuronFile(filename, header, headerSize, triggers);

% save out all of the data
for i=1:length(spikes)
  neuronsFile.addNeuron(electrodes(i), cell_ids(i), spikes{i}, length(spikes{i}));
end

% close it up
neuronsFile.close;

% update user
disp('done.');