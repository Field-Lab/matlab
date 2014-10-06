function p = Plot(d, nElectrodes)
%PLOT Plot class constructor.
%   p = PLOT(Data d) creates a plot object from Data object d,
%   where d contains spike vectors and projections.
%   present on a single electrode.
%
%   p = PLOT(Plot d) returns a copy of Plot object d.
%
%   tamachado@salk.edu 1/23/08

if nargin < 2
    nElectrodes = 512;
end


% Copy constructor
if isa(d,'Plot')
   p = d;

% Construct using a Data object
elseif isa(d, 'Data')
   p.data   = d;
   p.figure = cell(nElectrodes);
   p = class(p,'Plot');
   
% Otherwise generate an error
else
   disp('Warning: Invalid Constructor Call, No Object Created!')
end