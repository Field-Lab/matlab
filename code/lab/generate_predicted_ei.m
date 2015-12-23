function [axon,soma,axon_dists,soma_dists] = generate_predicted_ei(axon,positions, varargin)
% generate_predicted_ei     read in an axon and predict the EI
%
% usage:  ei = generate_predicted_ei(axon,positions, varargin)
%
% arguments:     axon - Nx2 matrix of axon path, starting at soma
%           positions - positions of electrodes
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     ei - E-length vector of EI intensity on each electorde
%
%
% optional params, their default values, and what they specify:
%
% method          	'basic'       	how to compare
%                                       'basic' - ???
%
%
% 2010-01  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('method', 'basic');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


switch params.method
    case 'basic'
        
        axon_dists = zeros(size(positions,1),1);
        soma_dists = zeros(size(positions,1),1);

        for ee =1:size(positions,1)
            x = positions(ee,1);
            y = positions(ee,2);
            % compute distance from each point to the axon
            axon_dists(ee) = sqrt(min(   (x-axon(:,1)).^2 + (y-axon(:,2)).^2   ));
            soma_dists(ee) = sqrt((x-axon(1,1)).^2 + (y-axon(1,2)).^2);
        end

        % make axon
        axon = (1./(axon_dists + 10)).^1;
        %axon = (axon_dists + 15).^-0.5;
        %axon = (axon_dists + 10).^-1.8;
        
        % make soma
        soma = 1./(5 + soma_dists.^2 );
end
