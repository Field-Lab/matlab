function maps = load_stimpixmaps(piece)

datadir = [server_data_path piece '/'];
datacontents = dir([datadir piece '_*']);

maps = struct([]);
for i = 1:length(datacontents)
    if ~datacontents(i).isdir, continue; end
    thisdir = datacontents(i).name;
    
    thisdirmaps = dir([datadir thisdir '/map-*.txt']);
    if ~isempty(thisdirmaps)
        maps(end+1).name = thisdir;
        maps(end).maps = {};
        for j = 1:length(thisdirmaps)
            thismap = thisdirmaps(j).name;
            maps(end).maps{end+1} = load_stimpixmap([datadir thisdir '/' thismap]);
        end
    end
end