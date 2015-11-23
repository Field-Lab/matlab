close all
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2009-04-13-3/d08_12_16_19-20-norefit/data008-from-data008_data012_data016_data019_data020/data008-from-data008_data012_data016_data019_data020.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2009-04-13-3/d08_12_16_19-20-norefit/data008-from-data008_data012_data016_data019_data020/data008-from-data008_data012_data016_data019_data020.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2009-04-13-3/d08_12_16_19-20-norefit/data008-from-data008_data012_data016_data019_data020/data008-from-data008_data012_data016_data019_data020.sta';
% datarun.names.rrs_neurons_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.neurons';
% datarun.names.rrs_params_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.params';
% datarun.names.rrs_sta_path = '/Volumes/Analysis/2015-10-06-8/data001/data001.sta';

vision_id = 7339;
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
datarun.names.rrs_ei_path = '/Volumes/Analysis/2009-04-13-3/d08_12_16_19-20-norefit/data008-from-data008_data012_data016_data019_data020/data008-from-data008_data012_data016_data019_data020.ei';
datarun=load_ei(datarun, vision_id);
% plot_ei_scroll(datarun,7339)


positions = datarun.ei.position;
[xc,yc] = getElectrodeCoords512(); 
 
number_of_frames = size(datarun.ei.eis{get_cell_indices(datarun, vision_id)},2);
for i = 1:number_of_frames
    figure
scatter(xc,yc, 4*ones(size(xc)), 'k','filled')
hold on
plot_ei(datarun,vision_id, 'frame_number',i,'foa', gca)
axis off
axis equal
F(i) = getframe(gcf);
close(gcf)
end

% movie2avi(F, 'EImovie.avi')


v = VideoWriter('EI movie');
v.FrameRate = 30%number_of_frames/3; 
open(v)
writeVideo(v,F);
close(v)
