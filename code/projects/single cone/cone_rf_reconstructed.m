function cone_rf = cone_rf_reconstructed(rf_size,cone_weights,cone_centers, varargin)
% cone_rf_reconstructed     image of the cone sampling of a single RF based on extracted cone weights
%
% usage:  cone_rf = cone_rf_reconstructed(rf_size,cone_weights,cone_centers, varargin)
%
% arguments:
%             rf_size - 2-length vector of height and width of an RF
%        cone_weights - N-length vector of cone weights
%        cone_centers - Nx2 matrix, x,y coordinates of each cone center point
%            varargin - struct or list of optional parameters (see below)
%
% outputs:
%       if cone_types == [],
%           cone_rf - YxX matrix of the cone sampling of the specified cell
%
%       otherwise,
%           cone_rf - YxXx3 matrix of the cone sampling of the specified cell, 
%                       with a different cone type in each color channel
%
%
% optional params, their default values, and what they specify:
%
% cone_types      	[]          optionally color cones according to their type
% cones             1:length(cone_weights,1)
%                               which cones to use.  N-length vector of cone IDs
%
%
% gauthier   2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('cone_types', []);
p.addParamValue('cones', 1:length(cone_weights));


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% BODY OF THE FUNCTION

% ensure params.cones is a row vector
if size(params.cones,2) == 1
    params.cones = params.cones';
end
    
% initialize with zeros
cone_rf = zeros(rf_size);

if isempty(params.cone_types)

    % go through each cone
    for cc = params.cones
        % enter weight in the approximately correct pixel
        cone_rf( min(max( round(cone_centers(cc,2)),1),rf_size(1)) ,...
            min(max( round(cone_centers(cc,1)),1),rf_size(2)) ) = cone_weights(cc);
    end

else

    % go through each cone
    for cc = params.cones

        % note which colors to enter
        switch params.cone_types(cc)
            case 'L'; cols = 1;
            case 'M'; cols = 2;
            case 'S'; cols = 3;
            case 'U'; cols = 1:3;
        end

        % enter weight in the approximately correct pixel
        cone_rf(round(cone_centers(cc,2)),round(cone_centers(cc,1)),cols) = cone_weights(cc);
    end

end

