function save_cones(path2save, mapName)

global datarun stim cellType myCells

height = datarun.stimulus.stixel_height * datarun.stimulus.field_height;
width = datarun.stimulus.stixel_width * datarun.stimulus.field_width;

myMap=zeros(height,width);


for i=1:length(stim.coord)
    myMap(stim.coord(i,1) - stim.stimarea:stim.coord(i,1) + stim.stimarea, ...
        stim.coord(i,2) - stim.stimarea:stim.coord(i,2) + stim.stimarea) = i;
end

figure
imagesc(myMap);

if ~exist(path2save,'dir')
    mkdir(path2save)
end

clear info
info.name = datarun.names.rrs_prefix;
info.cellType = cellType;
info.cells = myCells;

dlmwrite([path2save mapName, '.txt'], myMap, 'delimiter', '\t', 'newline', 'pc');
save([path2save mapName '_info'],'stim', 'info')



