function plotBundle(positions, arrayVoltages, path)
f = figure; set(f,'Position',[100 360 1000 550]);
set(f,'Color','white');
cla;
scatter(positions(:,1),positions(:,2),350,arrayVoltages,'filled');
axis off; axis image; colorbar;
caxis([-50 -10]);

hold on; scatter(positions(path,1),positions(path,2),350,[0.5 0.5 0.5], 'filled');
%         text(positions(stimChan,1),positions(stimChan,2),'stimulating electrode')
pause(0.001);
end