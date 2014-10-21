clear;

%% Parameters
stimulatedNeurons = [4726 4923 4951 5012 5221 5237 5327 5401 5432 5551 5671 5793 ...
                     5821 5866 5882 6031 6079 6271 6334 6361 6410 6692 ...
                     6767 6887 541 571 692 783 976 1098 1187 1201 1564 1696 ...
                     1966 2131 2193 2296 2356 2401 2461 2657 2746 4171 4426 4441 ...
                     4561 4698 4726 5237 5327];
             
%% Plotting the result

figuresFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data002/figures/spikes_vs_dist2d';
statsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data002/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data002/vision_processing/data002';
% 2013-04-04-0, data002
spotCenterPosition = [660;300;225;210];
for kk=1:10
    plotElicitedSpikesVsDistance2Px(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...
        'stimToPlot',kk,...
        'figNumber',2,...
        'stimulatedNeuronsList',stimulatedNeurons);
end

figuresFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/figures/spikes_vs_dist2d';
statsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/vision_processing/data004';
% 2013-04-04-0, data004
spotCenterPosition = [615;-130;55;95];
for kk=1:10
    plotElicitedSpikesVsDistance2Px(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...
        'stimToPlot',kk,...
        'figNumber',2,...
        'stimulatedNeuronsList',stimulatedNeurons);
end

figuresFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/figures/spikes_vs_dist2d';
statsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data004/vision_processing/data004';
% 2013-04-04-0, data004
spotCenterPosition = [615;-130;-255;-90];
for kk=1:10
    plotElicitedSpikesVsDistance2Px(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...
        'stimToPlot',kk,...
        'figNumber',2,...
        'stimulatedNeuronsList',stimulatedNeurons);
end

figuresFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data005/figures/spikes_vs_dist2d';
statsFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data005/statistics/activation';
dataFolder = '/media/MEA_PROCESSED_5/2013-04-04-0/data/data005/vision_processing/data005';
% 2013-04-04-0, data005
spotCenterPosition = [435;210;-255;-90];
for kk=1:10
    plotElicitedSpikesVsDistance2Px(dataFolder,statsFolder,figuresFolder,spotCenterPosition,...
        'stimToPlot',kk,...
        'figNumber',2,...
        'stimulatedNeuronsList',stimulatedNeurons);
end
