function save_all(datarun, cones, myMap, path2save, mapName, stim)



if ~exist(path2save,'dir')
    mkdir(path2save)
end

info.name = datarun.names.rrs_prefix;
info.cellType = cellType;
info.cells = myCells;

dlmwrite([path2save mapName, '.txt'], myMap, 'delimiter', '\t', 'newline', 'pc');
save([path2save mapName '_info'],'stim', 'info', 'cones')

close(hnew)

figure(tmp);

hSavedInfo=uicontrol('style','text', 'Units', 'Normalized','position',[0.32 0.08 0.32 0.1],...
    'string',{'saved to', path2save, ['file: ' mapName '.txt'], [mapName '_info.mat']},'fontsize',16, 'fontweight', 'bold','foregroundcolor','b');

