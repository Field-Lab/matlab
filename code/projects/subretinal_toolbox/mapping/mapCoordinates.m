function [arrayToMEA MEAToArray] = mapCoordinates(pixelCoordsMEA, pixelIndexArray, ...
                                       deviceType, flipped)
% [ARRAYTOMEA MEATOARRAY] = mapCoordinates(pixelCoordsOnMEA, ...
%                                             pixelIndexArray, deviceType, flipped)
%
% This function computes the rotation and translation matrix neccessary to
% go from MEA-based coordinates to pixel coordinates, and vice-versa from
% pixel coordinates to MEA-based coordinates. 
% The mapping is computed using Horn's absolute orientation method.
%
% Parameters:
%   - pixelCoordsOnMea: a n*2 matrix of coordinates of array pixels on the
%   MEA, specified in microns. The more pixels, the more accurate the
%   mapping will be.
%   - pixelIndexArray: a vector with n elements. Each element is the number
%   of the pixel whose position on the MEA is pixelCoordsOnMea(kk,:). Pixel
%   number can be found in the array_lattice files for each device type.
%   - deviceType: string defining what the type of the device is. It can
%   either be 'small' or 'medium'.
%   - flipped: a flag which is set to true if the device should be flipped.
%
% Optional parameters:
%   []
%
% Returns:
%   - arrayToMEA: a struct with fields arrayToMEA.R and arrayToMEA.t, where
%   the first one is the rotation matrix to go from array coordinates to
%   MEA coordinates, and the second one is the translation vector.
%   - MEAToArray: similar structure to go from MEA to array coordinates. 

% Version: 1.0 - 05/31/2013
% Author: Georges Goetz - Stanford University

%% Load the pixel position information

if ~exist('flipped','var')
    flipped = false;
end

switch deviceType
    case 'medium'
        if flipped
            pixelCoordsArray = dlmread('med_array_lattice_flipped-lr.txt');
        else
            pixelCoordsArray = dlmread('med_array_lattice.txt');
        end
    case 'small'
        if flipped
            pixelCoordsArray = dlmread('small_array_lattice_flipped-lr.txt');
        else
            pixelCoordsArray = dlmread('small_array_lattice.txt');
        end
    case 'large'
        error('Not yet implemented');
end
pixelCoordsArray = pixelCoordsArray(pixelIndexArray,:);

%% Compute the mapping

% From array to MEA
regParams = absor(pixelCoordsArray.', pixelCoordsMEA.');
arrayToMEA.R = regParams.R;
arrayToMEA.t = regParams.t;

% From MEA to array
regParams = absor(pixelCoordsMEA.', pixelCoordsArray.');
MEAToArray.R = regParams.R;
MEAToArray.t = regParams.t;

end % mapMEAArrayCoordinates