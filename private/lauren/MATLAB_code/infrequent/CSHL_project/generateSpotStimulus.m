function stimulus = generateSpotStimulus(spotCenter, spotRadius, stimSize)

% generates matrix of grayscale values with spot centered at spotCenter,
% with spotRadius (units = microns) with overall stimulus size stimSize
% 
% stimulus: 
%    first index = increasing x position (1:stimSize(1))
%    second index = increasing y positioon (1:stimSize(2))
%% for testing
% clear all
%  spotCenter = [50 100];
%  spotRadius = 10;
%  stimSize = [300 200]; %width x height

%%
th = 0:0.01:2*pi;

xCircle = cos(th)*spotRadius + spotCenter(1);
yCircle = sin(th)*spotRadius + spotCenter(2);

%plot(xCircle, yCircle)
%axis equal

xGrid = 1:stimSize(1);
yGrid = 1:stimSize(2);
[xGrid, yGrid] = meshgrid(xGrid, yGrid); %xGrid = repeated row vectors, yGrid = repeated column vectors

%coordinates of stimulus: 
%first index (moving down through rows) = increasing y values,
%second index (moving over through columns) = increasing x values
stimulus = inpolygon(xGrid, yGrid, xCircle, yCircle);

%flip so that first index = increasing x position, second index =
%increasing y position

stimulus = double(stimulus');
% figure
% imagesc(stimulus);
% colormap(gray)
% axis equal

