clear;

%% Parameters

% 2013-05-22-0
deviceType = 'medium';
pixelIndexArray = 1:37; % Index of the pixels that were stimulated in the array
stimulatedNeurons = [437 481 678 916 1022 1051 1262 1456 1743 3946 4141 4142 ...
    4216 4292 4322 4367 4369 4427 4442 4653 4666 4682 4727 4805 4951 4997 ...
    5177 5311 5359 5416];
figuresFolder = '/media/MEA_PROCESSED_5/2013-05-22-0/data/data002/figures/spikes_vs_dist';
statsFolder = '/media/MEA_PROCESSED_5/2013-05-22-0/data/data002/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-05-22-0/data/data002/vision_processing/data002';

%% Getting the spot center positions

% Pixels used for the mapping of coordinates
pixelIndex = [1;4;16;34;37];
pixelCoordsOnMEA = [-315 -150; -405 270; 105 -270; 435 -10; 345 450];
[arrayToMEA MEAToArray] = mapCoordinates(pixelCoordsOnMEA, pixelIndex, deviceType);

% Loading the right pixel coordinates file and use the mapping to position
% on the array.
switch deviceType
    case 'medium'
        spotCenterPosition = dlmread('med_array_lattice.txt');
    case 'small'
        spotCenterPosition = dlmread('small_array_lattice.txt');
    case 'large'
        error('Not yet implemented');
end
spotCenterPosition = spotCenterPosition(pixelIndexArray,:);
nSpots = length(spotCenterPosition);
for kk=1:nSpots
    spotCenterPosition(kk,:) = (arrayToMEA.R*spotCenterPosition(kk,:).' + arrayToMEA.t).';
end
             
%% Plotting the result

p = zeros(nSpots,3);

for kk=1:nSpots
    p(kk,:) = plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition(kk,:),...
        'stimToPlot',kk,...
        'figNumber',2,...
        'stimulatedNeuronsList',stimulatedNeurons);
end

%% Plotting all the fits together, getting all the 50% and 10% thresholds

f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);
f_inv = @(p,y) (-p(2)^2.*log(y./p(1))).^(1/p(3));
x = linspace(0,1000,5000);
thresholds = zeros(nSpots,2);

for kk=1:nSpots
    y = f(p(kk,:),x);
    thresholds(kk,1) = x(find(y<y(1)*0.5,1,'first'));
    thresholds(kk,2) = x(find(y<y(1)*0.1,1,'first'));
end

fh = figure(2); clf; set(fh,'color','white');
hold on;

y_ml = logspace(-3,0,5000);
mean_line = zeros(size(x));
for kk=1:nSpots
    mean_line = mean_line + f_inv([1 p(kk,2:3)],y_ml);
end
mean_line = mean_line/nSpots;
line(mean_line,y_ml,'color','r','linewidth',2); % for the legend

for kk=1:nSpots
    line(x,f([1 p(kk,2:3)],x),'color','k');
    mean_line = mean_line + f([1 p(kk,2:3)],x);
end
line(mean_line,y_ml,'color','r','linewidth',2); % plot again on top

title('All fits, NIR stimulation, medium 3-d devices with connected returns','fontsize',12)
xlabel('Distance to the center of the spot');
ylabel('Normalized number of spikes per trial')
axis([0 1000 0 1])
legend('Average fit')
