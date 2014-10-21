function [p, f] = predictResponseToLightSpot(dataFolder, figuresFolder, lightSpotCenterPos, lightSpotDiameter, varargin)
% predictResponseToLightSpot(dataFolder, lightSpotCenterPos, lightSpotDiameter)
%
% This function predicts what the response of a population of neurons whose
% receptive fields have been calculated to a light spot will be.
% The position of the receptive fields is estimated from the position of
% the EIs, doing a least squares fit of the receptive fields position onto
% the EI positions.
%
% Parameters:
%   - dataFolder: folder in which a .neurons and a .params files can be
%   found
%   - figuresFolder: folder in which the plots should be saved.
%   - lightSpotCenterPos: position of the center of the light spot in
%   microns. The point (0,0) corresponds to the center of the
%   microelectrode array.
%   - lightSpotDiameter: diameter of the light spot, in microns.
%
% Optional 'parameter'/value pairs:
%   TODO
%
% Returns:
%   - p: parameters of the function that was fit onto the datapoints.
%   - f: the function that was fit onto the datapoints. Typically, 
%   f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);

% Version: 0.0.2 - 4/02/2013
% 

%% Default values

if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end

%% Reading the input arguments

savePlot = true;
titleString = sprintf('Predicted response of RGCs to a %dum diameter spot of light', lightSpotDiameter);
figNumber = 1;
imageFormat = 'png';
neuronList = [];
stixelToUm = 134;
                 
% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'titlestring'
            titleString = varargin{kk*2};
        case 'fignumber'
            figNumber = varargin{kk*2};
        case 'imageformat'
            imageFormat = varargin{kk*2};
        case 'saveplot'
            savePlot = varargin{kk*2};
        case 'neuronlist'
            neuronList = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Making sure folders ends by '\' or '/'
if dataFolder(end:end)~=filesep
    dataFolder = [dataFolder filesep];
end
if figuresFolder(end:end)~=filesep
    figuresFolder = [figuresFolder filesep];
end
if size(lightSpotCenterPos,2) > size(lightSpotCenterPos,1)
    lightSpotCenterPos = lightSpotCenterPos.';
end

if ~exist(figuresFolder,'dir')
    mkdir(figuresFolder);
end

%% Linking to the files

if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

% Finding the neuron file and the param file
contentsDataFolder = dir(dataFolder);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isParamFile = strfind(contentsDataFolder(kk).name,'.params');
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuron_path =  [dataFolder contentsDataFolder(kk).name];
            end
        end
    end
    if isParamFile
        param_path = [dataFolder contentsDataFolder(kk).name];
    end
end
% Linking to the files
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
paramFile=edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

%% Reading the information for each neuron

neuronInfo = double(neuronFile.getIDList());
if ~isempty(neuronList)
    neuronInfo = intersect(neuronList, neuronInfo);
end
if size(neuronInfo,1)<size(neuronInfo,2)
    neuronInfo = neuronInfo.';
end
nNeurons = length(neuronInfo);

% Storing in order: neuronID, xPos center rec. field, yPos, xSigma, ySigma,
% ThetaRF, EIx0, EIy0
neuronInfo = [neuronInfo repmat(zeros(size(neuronInfo,1),1),1,6)];
for kk=1:nNeurons
    neuronInfo(kk,2) = paramFile.getDoubleCell(neuronInfo(kk,1),'x0');
    neuronInfo(kk,3) = paramFile.getDoubleCell(neuronInfo(kk,1),'y0');
    neuronInfo(kk,4) = paramFile.getDoubleCell(neuronInfo(kk,1),'SigmaX');
    neuronInfo(kk,5) = paramFile.getDoubleCell(neuronInfo(kk,1),'SigmaY');
    neuronInfo(kk,6) = paramFile.getDoubleCell(neuronInfo(kk,1),'Theta');
    neuronInfo(kk,7) = paramFile.getDoubleCell(neuronInfo(kk,1),'EIx0');
    neuronInfo(kk,8) = paramFile.getDoubleCell(neuronInfo(kk,1),'EIy0');
end

% Removing the bad neurons from the list
badNeurons = false(nNeurons,1);
for kk=1:size(neuronInfo,2)
    badNeurons = badNeurons | logical(isnan(neuronInfo(:,kk)));
end
neuronInfo(badNeurons,:) = [];
nNeurons = length(neuronInfo);

% Then, estimate the true receptive field position for each of the neurons
% For how it is done, look up Absolute Orientation, Horn's method
% The neuron list needs to be clean for this method to work, it's quite
% sensitive to outliers and will fail if there are too many of them (so
% clean up STAs!).
% Also, the columns of the data structure are magic numbers, which is
% sorta evil.
[regParams,Bfit,ErrorStats]=absor(neuronInfo(:,2:3).',neuronInfo(:,7:8).','doScale',1);

% Finally, adjust the receptive fields positions and sizes
neuronInfo(:,2:3) = (regParams.s*regParams.R*neuronInfo(:,2:3).' + repmat(regParams.t,1,length(neuronInfo))).';
neuronInfo(:,4:5) = (regParams.s)*neuronInfo(:,4:5);

%% Compute the overlap and distance between receptive field and spot for each neuron

overlapRF = zeros(nNeurons,1);
distToCenter = zeros(nNeurons,1);
for kk=1:nNeurons
    RFCenter= neuronInfo(kk,2:3);
    RFShape = neuronInfo(kk,4:5);
    RFAngle = neuronInfo(kk,6);
    overlapRF(kk) = computeOverlapReceptiveFieldLightSpot(lightSpotCenterPos, ...
        lightSpotDiameter, RFCenter, RFShape, RFAngle);
    distToCenter(kk) = computeDistanceRFToSpotCenter(lightSpotCenterPos, RFCenter);
end

%% Fit a (quasi) Gaussian to the result

% Fit a (quasi) Gaussian onto the datapoints
f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);

% Specifying fitting options
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 10 1.5],'Upper',[1.2*max(overlapRF) 1000 2]);
st_ = [10 200 1.5];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*exp(-abs(x^c)/b^2)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b', 'c'});

% Doing the fitting
cf_ = fit(distToCenter,overlapRF,ft_,fo_);
p = coeffvalues(cf_);

%% Plot the result

fh = figure(figNumber); clf; set(fh,'color','white');
hold on
scatter(distToCenter, overlapRF);
line(linspace(0,max(distToCenter),400),f(p,linspace(0,max(distToCenter),400)),...
    'color','r')

title(titleString,'fontsize',12);
xlabel('Distance to the center of the spot, um');
ylabel('Response of the RGCs, a.u.')
grid on
axis([0 min(max(500,4*lightSpotDiameter),900) 0 1.1*max(max(overlapRF),f(p,0))]);

%% Saving

if savePlot
    saveas(fh,sprintf('%sresponse_spotCenter%d-%d_spotDiam%d',figuresFolder,...
        lightSpotCenterPos(1),lightSpotCenterPos(2),lightSpotDiameter),imageFormat);
end

end % predictResponseToLightSpot

function overlap = computeOverlapReceptiveFieldLightSpot(lightSpotCenterPos, ...
        lightSpotDiameter, RFCenter, RFShape, RFAngle)
% This function computes the overlap between a spot of light of a given
% diameter and a receptive field modelled by a Gaussian distribution of
% know covariance. 
% The receptive field is scaled so that the value of the Gaussian at the
% center is always 1 (so the area under the Gaussian is not 1).
% The result is returned in arbitrary units.
    
% Creating estimation grid
nPointsEstimation = round(sqrt(10000));
xMin = lightSpotCenterPos(1) - lightSpotDiameter;
xMax = lightSpotCenterPos(1) + lightSpotDiameter;
yMin = lightSpotCenterPos(2) - lightSpotDiameter;
yMax = lightSpotCenterPos(2) + lightSpotDiameter;
[XX, YY] = meshgrid(xMin:((xMax-xMin)/nPointsEstimation):xMax, ...
    yMin:((yMax-yMin)/nPointsEstimation):yMax);

% Creating the spot
spotValue = double(sqrt((XX - lightSpotCenterPos(1)).^2 + (YY - lightSpotCenterPos(2)).^2)...
    <lightSpotDiameter/2);
spotValue = spotValue/sum(spotValue(:));

% Creating the RF
% display(sprintf('SigmaX: %d, SigmaY:%d',round(RFShape(1)),round(RFShape(2))));
RotMatrix = [cos(RFAngle) -sin(RFAngle); sin(RFAngle) cos(RFAngle)];
RFCovMatrix = RotMatrix*diag(RFShape.*RFShape)*RotMatrix.';

RFValue = reshape(mvnpdf([XX(:) YY(:)], RFCenter, RFCovMatrix),size(XX,1),size(XX,1));
scalingFactor = 1/mvnpdf(RFCenter, RFCenter, RFCovMatrix);
RFValue = RFValue*scalingFactor;

% Estimate overlap
overlap = sum(sum(spotValue.*RFValue));

end % computeOverlapReceptiveFieldLightSpot

function distToCenter = computeDistanceRFToSpotCenter(lightSpotCenterPos, RFCenter)
% This function computes the distance between a receptive field center and
% the center of a light spot.

if size(lightSpotCenterPos,1)<size(lightSpotCenterPos,2)
    lightSpotCenterPos = lightSpotCenterPos.';
end
if size(RFCenter,1)<size(RFCenter,2)
    RFCenter = RFCenter.';
end

distToCenter = norm(lightSpotCenterPos - RFCenter);

end % computeDistanceRFToSpotCenter