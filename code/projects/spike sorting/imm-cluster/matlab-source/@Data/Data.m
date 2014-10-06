function d = Data(path, params)
%DATA Data class constructor.
%   d = DATA(path, params) open up projections file at path. Subsequent
%   calls to getStruct(electrode) return a struct containing
%   data from a specific electrode.
%
%
%   tamachado@salk.edu 1/23/08



% Copy constructor
if isa(path,'Data')
   d = path;
   
% Open file stream to read from prj file
elseif isa(path, 'char')
   prjObject = edu.ucsc.neurobiology.vision.matlab.ReadProjections(path);
   d.prjObject     = prjObject;
   d.prjStruct     = [];
   d.badElectrodes = false(1, params.nElectrodes);
   d = class(d,'Data');

% Otherwise generate an error
else
   disp('Warning: Invalid Constructor Call, No Object Created!')
end