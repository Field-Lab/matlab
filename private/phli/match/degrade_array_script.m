datarun = load_data('2011-07-05-4/data003');
datarun = load_neurons(datarun);
datarun = load_ei(datarun, []);


%%
positions = datarun.ei.position;
arraytype = size(positions, 1);


%% Determine inline direction
start = 4;
neighbors = get_ei_neighbors(start, arraytype);
nstruct = get_lazyhex_ei_neighbors(start, neighbors, positions);

inline = setdiff(nstruct.inline, start);
inlinediff = positions(start,:) - positions(inline(1),:);


%% Get seed electrodes for every other line
deglineseeds = find(positions(:,2) == positions(start,2));
[~, ind] = sort(positions(deglineseeds,1));
deglineseeds = deglineseeds(ind);


%% Get each line of electrodes, sort according to position in line, degrade
deglines = cell(length(deglineseeds),1);
for i = 1:length(lines)
    seed = deglineseeds(i);
    
    % Have to be a little fuzzy for some reason
    line = find(positions(:,1) > positions(seed,1)-1 & positions(:,1) < positions(seed,1)+1);
    
    % Sort
    [~, ind] = sort(positions(line,2));
    line = line(ind);
    
    % Degrade and save
    degline = line(1:2:end);
    deglines{i} = degline;
end


%%
degelecs = cell2mat(deglines);
degpos = positions(degelecs,:);
plot(degpos(:,1), degpos(:,2), '.')


%%
save ../matlab-standard/private/phli/match/degelecs519 degelecs