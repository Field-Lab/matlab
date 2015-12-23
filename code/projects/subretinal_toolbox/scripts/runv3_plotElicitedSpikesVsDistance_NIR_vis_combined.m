clear;

%% Plotting the #spikes vs dist

p = zeros(4,3);

% 2013-03-13-0, data003 parameters: NIR spots
stimulatedNeurons = []; 
figuresFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data003/figures/spikes_vs_dist';
statsFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data003/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data003/vision_processing/data003';
spotCenterPosition = [105 -390;
                      -165 150].';
trials_NIR = [1 13];

for kk=1:length(trials_NIR)
    p(kk,:) = plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition(:,kk),...
        'stimToPlot',trials_NIR(kk),...
        'figNumber',1,...
        'stimulatedNeuronsList',stimulatedNeurons,...
        'savePlot',false);
end

% 2013-03-13-0, data006 parameters: visible spots
stimulatedNeurons = []; 
figuresFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data006/figures/spikes_vs_dist';
statsFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data006/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-03-13-0/data/data006/vision_processing/data006';
spotCenterPosition = [105 -390;
                      -165 150].';
trials_vis = [1 9];

for kk=1:length(trials_vis)
    p(kk+length(trials_NIR),:) = plotElicitedSpikesVsDistance(dataFolder,statsFolder,figuresFolder,spotCenterPosition(:,kk),...
        'stimToPlot',trials_vis(kk),...
        'figNumber',1,...
        'stimulatedNeuronsList',stimulatedNeurons,...
        'savePlot',false);
end

%% Plotting

npoints = 5000;
maxdist = 2000;
nLines = size(p,1);
colorList = {'r';'b';'--r';'--b'};

f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);
f_inv = @(p,y) (-p(2)^2.*log(y./p(1))).^(1/p(3));
x = linspace(0,maxdist,npoints);

fh = figure(2); clf; set(fh,'color','white');
hold on;

for kk=1:nLines
    plot(x,f([1 p(kk,2:3)],x),colorList{kk},'linewidth',2);
end

title('Comparison of visible and NIR spread, 2 spots','fontsize',12)
xlabel('Distance to the center of the spot');
ylabel('Normalized number of spikes per trial')
axis([0 1000 0 1])
legend('NIR spot 1','NIR spot 2','Visible spot 1','Visible spot 2')
xlim([0 800])