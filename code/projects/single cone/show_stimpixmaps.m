function show_stimpixmaps(mapsets)

for i = 1:length(mapsets)
    thismapset = mapsets(i);
    nummaps = length(thismapset.maps);
    
    for j = 1:nummaps
        sanesubplot(nummaps, length(mapsets), [j, i]);
        imagesc(thismapset.maps{j});
        axis equal tight
        set(gca, 'XTick', get(gca, 'XLim') - 0.5);
        set(gca, 'YTick', get(gca, 'YLim') - 0.5);
        
        if j == 1
            title(thismapset.name, 'Interpreter', 'none');
        end
    end
end