function [m, err] = ModelFile(path, min, oldModel)
%MODELFILE ModelFile class constructor.
%  m = MODELFILE(path, data, nBest) opens a model file at at path.
%  Subsequent calls to add() will add clusters to the model file.
%
%  Model files only save out clusters from the MAP
%
%   tamachado@salk.edu 4/3/08

% Copy constructor
if isa(path,'ModelFile')
   m = path;

% Open file stream to read from prj file
elseif isa(path, 'char')

   % Do this for each set of clusters we are saving out
   m.nObject   = edu.ucsc.neurobiology.vision.matlab.ReadModel(path, oldModel);
   
   err = m.nObject.getError();
   
   m.minSpikes = min;

   % Register this instance of the neurons file
   m = class(m,'ModelFile');

% Otherwise generate an error
else
   disp('Warning: Invalid Constructor Call, No Object Created!')
end