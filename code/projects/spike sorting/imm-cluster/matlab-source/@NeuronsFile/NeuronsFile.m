function n = NeuronsFile(path, data, min, nBest)
%NEURONSFILE NeuronFile class constructor.
%  n = NEURONSFILE(path, data, min, nBest) opens a neurons file at path.
%  Subsequent calls to add() will add clusters to the neurons file.
%
%
%   tamachado@salk.edu 1/31/08

% Copy constructor
if isa(path,'NeuronsFile')
   n = path;

% Open file stream to read from prj file
elseif isa(path, 'char') && isa(data, 'Data') 
   
   % Get the prj file object 
   prjObject = toJavaObject(data);
   prj = prjObject.getProjectionsFile;
   
   % Save the original path before appending an index to the end
   pO = path;
   
   for k = 1:nBest
       % Create/open the neurons file
       
       % If line is commented out, output only one file, otherwise output
       % many files
       if nBest > 1
           %path = [pO sprintf('%d', k)];
       end

       % Do this for each set of clusters we are saving out
       n.nObject{k}   = edu.ucsc.neurobiology.vision.matlab.ReadNeurons(path, prj);
       n.data{k}      = data;
       n.minSpikes{k} = min;
   end
   
   % Register this instance of the neurons file
   n = class(n,'NeuronsFile');

% Otherwise generate an error
else
   disp('Warning: Invalid Constructor Call, No Object Created!')
end