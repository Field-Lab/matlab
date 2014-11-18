%subplot moasics

javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

%   Initializing the data structure
%   The most important function is load data! load_data is expecting to be
%   directed to the file with the data000.bin, data000.neuron etc. These are
%   in the Analysis drive. The function is set up so that you can just type
%   'Date/data000'. This initializes the datarun structure. 
datarun=load_data('2007-03-27-1/data011-nwpca');

% Loading other information
%   Other information you might want including STAs, params, ei, neurons
%   (includes spike times) and polarities
datarun=load_sta(datarun);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
datarun=set_polarities(datarun);

% Let's look at the data! This is a classification run.

% Plot the mosaic for a cell type, in this case on parasols

figure;
subplot(4,2,1)
plot_rf_summaries(datarun,{1},'plot_fits',1,'label',0, 'coordinates','monitor')
title('ON Parasol Mosaic 2007-03-27-1')
axis equal
xlim([115 535])
ylim([125 410])

% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
set(gca,'YDir','reverse');

subplot(4,2,3)
plot_rf_summaries(datarun,{2},'plot_fits',1,'label',0, 'coordinates','monitor')
title('OFF Parasol Mosaic 2007-03-27-1')
axis equal
xlim([110 515])
ylim([130 410])

set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'YDir','reverse');

subplot(4,2,5)
plot_rf_summaries(datarun,{3},'plot_fits',1,'label',0, 'coordinates','monitor')
title('ON Midget Mosaic 2007-03-27-1')
axis equal
xlim([100 510])
ylim([130 355])

% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');

subplot(4,2,7)
plot_rf_summaries(datarun,{4},'plot_fits',1,'label',0, 'coordinates','monitor')
title('OFF Midget Mosaic 2007-03-27-1')
axis equal
xlim([110 485])
ylim([130 350])

% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');


datarun=load_data('2007-08-24-4/data001-nwpca');

% Loading other information
%   Other information you might want including STAs, params, ei, neurons
%   (includes spike times) and polarities
datarun=load_sta(datarun);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
datarun=set_polarities(datarun);

% Let's look at the data! This is a classification run.

% Plot the mosaic for a cell type, in this case on parasols
subplot(4,2,2)

plot_rf_summaries(datarun,{1},'plot_fits',1,'label',0, 'coordinates','monitor')
title('ON Parasol Mosaic 2007-08-24-4')
axis equal
xlim([125 510])
ylim([75.1 325])

% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');

subplot(4,2,4)
plot_rf_summaries(datarun,{2},'plot_fits',1,'label',0, 'coordinates','monitor')
title('OFF Parasol Mosaic 2007-08-24-4')
axis equal
xlim([130 490])
ylim([60 320])

set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'YDir','reverse');

subplot(4,2,6)
plot_rf_summaries(datarun,{3},'plot_fits',1,'label',0, 'coordinates','monitor')
title('ON Midget Mosaic 2007-08-24-4')
axis equal
xlim([150 495])
ylim([120 325])

set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'YDir','reverse');

set(gcf,'color','w');
