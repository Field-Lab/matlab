
figure('Color','w');
imagesc(0.1*mov_params.totalMaskAccept'  + cell_params.CellMasks{1}')
axis image
caxis([0,2]);