function imm = Imm(params, data, nElectrodes)
%IMM Imm class constructor.
%   imm = IMM(params, data) Do clustering using an infinite mixture model
%   This class fits PC projections stored in data to an infinite mixture
%   of gaussians model. The algorithm computes the best model to use by
%   sampling from a distribution of all solutions using MCMC techniques.
%   To cluster an electrode, use the clusterElectrode method.
%   The maximum a posteriori solution can be located using the getMAP
%   method. To extract the clusters from a solution found (e.g. the MAP)
%   use the getClusters method.
%
%
%   tamachado@salk.edu 1/28/08


%Example (after constructing an Imm object):
% The first step is to sample from the posterior using Imm.clusterElectrode
% Then we must fit the out of sample points and get the output using
% the Imm.getClusters method.
%
%[imm, mapIndex]  = clusterElectrode(imm, 1);
%cls  = getClusters(imm, mapIndex);


if nargin < 3
    nElectrodes = 512;
end

% Copy constructor
if isa(params,'Imm')
   imm = params;

% Set parameters for clustering and create a new Imm object
elseif isa(data, 'Data')
   imm.data     = data;                 %store projections/data
   imm.params   = params;               %store parameters for clustering
   imm.clusters = cell(nElectrodes);    %store generated custers
   imm.error    = zeros(1,nElectrodes); %record analyzed electrodes
   
   imm.SUCCESS_FIT     = 2;              %flag: success fitting out of sample pts          
   imm.SUCCESS_CLUSTER = 1;              %flag: success sampling from posterior
   
   imm = class(imm,'Imm');

% Otherwise generate an error
else
   disp('Warning: Invalid Constructor Call, No Object Created!')
end