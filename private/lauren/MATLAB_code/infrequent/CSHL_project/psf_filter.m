function filtered_stimulus = psf_filter(stimulus, cone_type)

% convolves estimated point spread function for ~12 mm eccentricity, based
% on data for 1 degree eccentricity provided by David Brainard
% arguments
%
%   stimulus: matrix of grayscale values (must be sampled at micron^2)
%   cone_type: 'L' or 'M' --> used to determine which psf to use
%
% output
%   filtered_stimulus: matrix of grayscale values (same dimensions as
%   stimulus)


%% for testing

% spotCenter = [250, 200];
% spotRadius = 30;
% stimSize = [500 500];
% 
% stimulus = generateSpotStimulus(spotCenter, spotRadius, stimSize);
% 
% cone_type = 'L';


%%
%in case stimulus is logical
stimulus = double(stimulus);

% apply point-spread function to stimulus
% 0.117 arc minutes per pixel, measured at 1 degree eccentricity (arc
% minute = 1/60 degree) --> 0.00195 degrees/pixel
if strcmpi(cone_type, 'L')
    psf = load('ap3lpsfSmall.txt');
elseif strcmpi(cone_type, 'M')
    psf = load('ap3mpsfSmall.txt');
end
%spsf = load('ap3spsfSmall.txt');

%scale to more peripheral eccentricity (roughly 12 mm vs. 0.27 mm)
%linear scale factor ~= 3.7x
%0.00195 degrees/pixel = 0.007215 degrees/pixel = 1.8038 microns/pixel
%(~250 microns/pixel)

%2D interpolation to get filters into correct dimensions (micron^2)
xOld = 0:1.8038:1.8038*(size(psf,1)-1);
yOld = 0:1.8038:1.8038*(size(psf,2)-1);
xNew = 0:1:1.8038*(size(psf,1)-1);
yNew = 0:1:1.8038*(size(psf,2)-1);
[xOldMesh yOldMesh] = meshgrid(xOld, yOld);
[xNewMesh yNewMesh] = meshgrid(xNew, yNew);

psf_interp = interp2(xOldMesh, yOldMesh, psf, xNewMesh, yNewMesh);

filtered_stimulus = conv2(stimulus, psf_interp, 'same');
%imagesc(filtered_stimulus)





