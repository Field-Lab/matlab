%subplot moasics

javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

%   Initializing the data structure
%   The most important function is load data! load_data is expecting to be
%   directed to the file with the data000.bin, data000.neuron etc. These are
%   in the Analysis drive. The function is set up so that you can just type
%   'Date/data000'. This initializes the datarun structure. 
datarun=load_data('2010-09-24-0/data001-nwpca');

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
ha = tight_subplot(1,2,[.05 .05],[.1 .1],[.1 .1]);
axes(ha(1));
plot_rf_summaries(datarun,{1},'plot_fits',1,'label',0, 'coordinates','monitor')
title('ON Parasol Mosaic')
axis equal
xlim([150 550])
ylim([0 500])

axis off
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'YDir','reverse');
hold on
% rectangle('Position', [110,130,395, 220])
axes(ha(2));
plot_rf_summaries(datarun,{2},'plot_fits',1,'label',0, 'coordinates','monitor')
title('OFF Parasol Mosaic')
axis equal
xlim([150 550])
ylim([0 500])

axis off

suptitle('2010-09-24-0/data001-nwpca')


set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'YDir','reverse');
% hold on
% rectangle('Position', [110,130,395, 220])
% axes(ha(3));
% plot_rf_summaries(datarun,{3},'plot_fits',1,'label',0, 'coordinates','monitor')
% title('ON Midget Mosaic')
% axis equal
% xlim([100 535])
% ylim([125 410])
% 
% axis off
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');
% hold on
% rectangle('Position', [110,130,395, 220])
% 
% axes(ha(4));
% plot_rf_summaries(datarun,{4},'plot_fits',1,'label',0, 'coordinates','monitor')
% title('OFF Midget Mosaic')
% axis equal
% xlim([100 535])
% ylim([125 410])
% 
% axis off
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');
% hold on
% rectangle('Position', [110,130,395, 220])

% datarun=load_data('2007-08-24-4/data001-nwpca');
% 
% % Loading other information
% %   Other information you might want including STAs, params, ei, neurons
% %   (includes spike times) and polarities
% datarun=load_sta(datarun);
% datarun=load_params(datarun);
% datarun=load_neurons(datarun);
% datarun=set_polarities(datarun);
% 
% % Let's look at the data! This is a classification run.
% 
% % Plot the mosaic for a cell type, in this case on parasols
% subplot(4,2,2)
% 
% plot_rf_summaries(datarun,{1},'plot_fits',1,'label',0, 'coordinates','monitor')
% title('ON Parasol Mosaic 2007-08-24-4')
% axis equal
% xlim([125 510])
% ylim([75.1 325])
% axis off
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');
% 
% subplot(4,2,4)
% plot_rf_summaries(datarun,{2},'plot_fits',1,'label',0, 'coordinates','monitor')
% title('OFF Parasol Mosaic 2007-08-24-4')
% axis equal
% xlim([125 510])
% ylim([75.1 325])
% axis off
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');
% 
% subplot(4,2,6)
% plot_rf_summaries(datarun,{3},'plot_fits',1,'label',0, 'coordinates','monitor')
% title('ON Midget Mosaic 2007-08-24-4')
% axis equal
% xlim([125 510])
% ylim([75.1 325])
% axis off
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% set(gca,'YDir','reverse');
% axis off
set(gcf,'color','w');
