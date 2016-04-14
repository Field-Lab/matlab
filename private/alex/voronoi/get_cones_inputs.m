
a = load('/Volumes/Analysis/2016-02-17-4/cone_data/manual/map_data001_manual_postexp_info.mat');
cones = round(a.cones/2);

mvi=load_movie('/Volumes/Analysis/stimuli/white-noise-xml/BW-2-6-0.48-11111-400x300-60.35.xml', datarun.triggers);
% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

filt_inputs(cnt,:) = filt_inputs(cnt,:) + (tmp_inputs - tmp_noise)*raw_sta(tmp_coord(2),tmp_coord(1));
cones_inputs = cell(1,length(cones));
tic
for j=0:duration-1
    if mod(j, 1000)==0
        j
    end
    F = round(mvi.getFrame(j).getBuffer);
    myFrames = reshape(F(1:3:end),width,height);    
    myFrames(myFrames==0) = -.48;
    myFrames(myFrames==1) = .48;
    for i = center_cones
        tmp = myFrames(cones(i,1)-1:cones(i,1)+1,cones(i,2)-1:cones(i,2)+1);
        cones_inputs{i} = [cones_inputs{i} tmp(:)];
    end
end
toc


save('/Users/alexth/Desktop/inputs_cell_541.mat', 'cones_inputs', 'refresh', 'duration', '-v7.3')
