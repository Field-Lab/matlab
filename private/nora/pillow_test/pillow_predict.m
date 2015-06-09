function xvalperformance = pillow_predict(fittedGLM, cid, testmovie_color, datarun, datarun_mas, trial_starts)

% Get Movie Params

stim_description = 'BW-20-1-0.48-11111';
disp(['Using ' stim_description ' XML file'])

% Master datarun to get RGB weights
master_idx         = find(datarun_mas.cell_ids == cid);
RGB = RGB_weights(datarun_mas,datarun_mas.cell_ids == cid);
master_idx         = find(datarun.cell_ids == cid);

% Turn RGB movie into greyscale movie
if strcmp(stim_description(1:2), 'BW')
    testmovie = squeeze(testmovie_color(:,:,1,:));
else
    testmovie=squeeze(RGB(1)*testmovie_color(:,:,1,:)+ ...
        RGB(2)*testmovie_color(:,:,2,:)+ ...
        RGB(3)*testmovie_color(:,:,3,:));
end

%% Get the test spikes

rawspikes = datarun.spikes{datarun.cell_ids == cid};
test_spikes = get_raster(rawspikes, trial_starts, 'plot', false);

%% Get neighbor spikes
paired_cells = fittedGLM.cellinfo.pairs;
for i = 1:length(paired_cells)
    neighbor_spikes{i} = get_raster(datarun.spikes{datarun.cell_ids == paired_cells(i)}, trial_starts, 'plot', false);
end

%%
xvalperformance = glm_predict(fittedGLM, testmovie, 'testspikes', test_spikes, 'neighborspikes', neighbor_spikes); 

end