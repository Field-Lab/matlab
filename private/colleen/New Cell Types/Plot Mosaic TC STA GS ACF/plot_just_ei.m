close all
clear

datarun.names.rrs_neurons_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.sta';
% datarun.names.rrs_neurons_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.neurons';
% datarun.names.rrs_params_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.params';
% datarun.names.rrs_sta_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.sta';

vision_id = [1130];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
datarun.names.rrs_ei_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.ei';

for j = 1:length(vision_id)
datarun=load_ei(datarun, vision_id(j), 'array_type', 519);
% plot_ei_scroll(datarun,7339)


positions = datarun.ei.position;
[xc,yc] = getElectrodeCoords519(); 
 
number_of_frames = size(datarun.ei.eis{get_cell_indices(datarun, vision_id(j))},2);
% for i = 1:number_of_frames
%     figure
% scatter(xc,yc, 4*ones(size(xc)), 'k','filled')
% hold on
% plot_ei(datarun,vision_id, 'frame_number',i,'foa', gca)
% axis off
% axis equal
% F(i) = getframe(gcf);
% close(gcf)
% end
% 
% % movie2avi(F, 'EImovie.avi')
% 
% 
% v = VideoWriter('EI movie');
% v.FrameRate = 30%number_of_frames/3; 
% open(v)
% writeVideo(v,F);
% close(v)

figure
scatter(xc,yc, 8*ones(size(xc)), 'k','filled')
hold on
plot_ei(datarun,vision_id(j), 'frame_number',0,'foa', gca)
axis off
axis equal
end

