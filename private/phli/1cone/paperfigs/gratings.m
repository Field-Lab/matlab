%% Basics
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};


%% 2011-12-13-2
piece = '2011-12-13-2';
conerun  = load_data(fullfile(piece, 'streamed', 'data000-0', 'data000-0'), staopts);
classrun = load_data(fullfile(piece, 'data001'), staopts);
gratingrun = load_data(fullfile(piece, 'data005'), loadopts);
gratingrun.mapconerun  = map_ei(conerun,  gratingrun);
gratingrun.mapclassrun = map_ei(classrun, gratingrun);
gratingrun = read_stim_lisp_output(gratingrun);
gratingrun.stimulus = parse_stim_rgbs(gratingrun.stimulus);

% Collect RGCs from both conerun and classification run
coneoffmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {4})}];
classoffmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {4})}];
offmidgets = union(coneoffmidg, classoffmidg);
coneoffpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {2})}];
classoffpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {2})}];
offparasols = union(coneoffpara, classoffpara);
coneonmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {3})}];
classonmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {3})}];
onmidgets = union(coneonmidg, classonmidg);
coneonpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {1})}];
classonpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {1})}];
onparasols = union(coneonpara, classonpara);

histx = 0:0.05:9.1;
for i = 1:length(offmidgets)
    offmidghists{i} = collect_spikehists(gratingrun, offmidgets(i), {'spatialperiods' 'urgbs'}, 'triggers', gratingrun.triggers(1:19:end));
    offmidgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offmidghists{i});
end
for i = 1:length(offparasols)
    offparahists{i} = collect_spikehists(gratingrun, offparasols(i), {'spatialperiods' 'urgbs'}, 'triggers', gratingrun.triggers(1:19:end));
    offparaf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offparahists{i});
end
for i = 1:length(onmidgets)
    onmidghists{i} = collect_spikehists(gratingrun, onmidgets(i), {'spatialperiods' 'urgbs'}, 'triggers', gratingrun.triggers(1:19:end));
    onmidgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), onmidghists{i});
end
for i = 1:length(onparasols)
    onparahists{i} = collect_spikehists(gratingrun, onparasols(i), {'spatialperiods' 'urgbs'}, 'triggers', gratingrun.triggers(1:19:end));
    onparaf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), onparahists{i});
end


%% 2013-03-11-1 data002 (can also analyze d05 d07 d17)
piece = '2013-03-11-1';
classrun = load_data(fullfile(piece, 'data000'), loadopts);
conerun = load_data(fullfile(piece, 'data001'), loadopts);
gratingrun = load_data(fullfile(piece, 'data002'), loadopts);
gratingrun.mapconerun  = map_ei(conerun,  gratingrun);
gratingrun.mapclassrun = map_ei(classrun, gratingrun);
gratingrun = read_stim_lisp_output(gratingrun);
gratingrun.stimulus = parse_stim_rgbs(gratingrun.stimulus);

% Collect RGCs from both conerun and classification run
coneoffmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {4})}];
classoffmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {4})}];
offmidgets = union(coneoffmidg, classoffmidg);
coneoffpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {2})}];
classoffpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {2})}];
offparasols = union(coneoffpara, classoffpara);
coneonmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {3})}];
classonmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {3})}];
onmidgets = union(coneonmidg, classonmidg);
coneonpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {1})}];
classonpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {1})}];
onparasols = union(coneonpara, classonpara);

histx = 0:0.05:8.1;
for i = 1:length(offmidgets)
    offmidghists{i} = collect_spikehists(gratingrun, offmidgets(i), {'spatialperiods' 'spatialphases' 'urgbs'}, 'triggers', gratingrun.triggers(1:17:end), 'histx', histx);
    offmidgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offmidghists{i});
end
offmidgoptf1s = cell2mat(cellfun(@(D)(squeeze(max(max(D)))), offmidgf1s, 'UniformOutput', false));
for i = 1:length(offparasols)
    offparahists{i} = collect_spikehists(gratingrun, offparasols(i), {'spatialperiods' 'spatialphases' 'urgbs'}, 'triggers', gratingrun.triggers(1:17:end), 'histx', histx);
    offparaf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offparahists{i});
end
offparaoptf1s = cell2mat(cellfun(@(D)(squeeze(max(max(D)))), offparaf1s, 'UniformOutput', false));


figure; hold on;
histedge = 0:0.5:6;
nmidg = histc(offmidgoptf1s(1,:), histedge);
npara = histc(offparaoptf1s(1,:), histedge);
hpara = bar(histedge, npara, 'histc');
hmidg = bar(histedge, nmidg, 'histc');
set(hmidg, 'FaceColor', 'k');
set(hpara, 'FaceColor', 'w');
axis([0 6.5 0 14]);


%% 2012-08-21-0 (probably could stand to clean up the classifications a bit)
piece = '2012-08-21-0';
% Not a real classrun, just another conerun; data000 is a real classrun,
% but I don't seem to have EIs at least not at home...
classrun = load_data(fullfile(piece, 'data001'), loadopts); 
conerun = load_data(fullfile(piece, 'data006'), loadopts);
gratingrun = load_data(fullfile(piece, 'data007'), loadopts);
gratingrun.mapconerun  = map_ei(conerun,  gratingrun);
gratingrun.mapclassrun = map_ei(classrun, gratingrun);
gratingrun = read_stim_lisp_output(gratingrun);
gratingrun.stimulus = parse_stim_rgbs(gratingrun.stimulus);

% Collect RGCs from both conerun and classification run
coneoffmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {4})}];
classoffmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {4})}];
offmidgets = union(coneoffmidg, classoffmidg);
coneoffpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {2})}];
classoffpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {2})}];
offparasols = union(coneoffpara, classoffpara);
coneonmidg  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {3})}];
classonmidg = [gratingrun.mapclassrun{get_cell_indices(classrun, {3})}];
onmidgets = union(coneonmidg, classonmidg);
coneonpara  = [gratingrun.mapconerun{ get_cell_indices(conerun,  {1})}];
classonpara = [gratingrun.mapclassrun{get_cell_indices(classrun, {1})}];
onparasols = union(coneonpara, classonpara);

histx = 0:0.05:8.1;
for i = 1:length(offmidgets)
    offmidghists{i} = collect_spikehists(gratingrun, offmidgets(i), {'spatialperiods' 'spatialphases' 'urgbs'}, 'triggers', gratingrun.triggers(1:17:end), 'histx', histx);
    offmidgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offmidghists{i});
end
offmidgoptf1s = cell2mat(cellfun(@(D)(squeeze(max(max(D)))), offmidgf1s, 'UniformOutput', false));
for i = 1:length(offparasols)
    offparahists{i} = collect_spikehists(gratingrun, offparasols(i), {'spatialperiods' 'spatialphases' 'urgbs'}, 'triggers', gratingrun.triggers(1:17:end), 'histx', histx);
    offparaf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), offparahists{i});
end
offparaoptf1s = cell2mat(cellfun(@(D)(squeeze(max(max(D)))), offparaf1s, 'UniformOutput', false));


figure; hold on;
histedge = 0:0.5:6;
nmidg = histc(offmidgoptf1s(1,:), histedge);
npara = histc(offparaoptf1s(1,:), histedge);
hpara = bar(histedge, npara, 'histc');
hmidg = bar(histedge, nmidg, 'histc');
set(hmidg, 'FaceColor', 'k');
set(hpara, 'FaceColor', 'w');
axis([0 6.5 0 14]);


%% 2012-09-13-2
piece = '2012-09-13-2';
conerun = load_data([piece '/streamed/data005/data005'], staopts);
conerun.offM = [1727 1816 1966 2026 2116 2386 2942 6991];
conerun.offP = [4441 5371 5551 6841];
conerun.rgcs = [conerun.offM conerun.offP];
% conerun.conepicks = [4 4 4 4 4 4 4 4 4 4 4 4];
gratingrun = load_data(fullfile(piece, 'data006'), loadopts);
gratingrun.map = map_ei(conerun, gratingrun);
gratingrun.offM = gratingrun.map(get_cell_indices(conerun, conerun.offM));
gratingrun.offP = gratingrun.map(get_cell_indices(conerun, conerun.offP));
gratingrun.rgcs = [gratingrun.offM gratingrun.offP];
gratingrun = read_stim_lisp_output(gratingrun);
gratingrun.stimulus = parse_stim_rgbs(gratingrun.stimulus);

histx = 0:0.05:8.1;

% Selected cells
midghists = collect_grating_spikehists(gratingrun, gratingrun.offM, gratingrun.triggers(1:17:end), 'histx', histx);
midghists = invertselect(midghists, @isempty);
for i = 1:length(midghists)
    midgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), midghists{i});
    midgoptf1s(:,i) = squeeze(max(max(midgf1s{i})));
end
parahists = collect_grating_spikehists(gratingrun, gratingrun.offP, gratingrun.triggers(1:17:end), 'histx', histx);
for i = 1:length(parahists)
    paraf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), parahists{i});
    paraoptf1s(:,i) = squeeze(max(max(paraf1s{i})));
end

% All OFF cells
midghists = collect_grating_spikehists(gratingrun, gratingrun.map(get_cell_indices(conerun, {4, 14, 16, 20, 24, 25})), gratingrun.triggers(1:17:end), 'histx', histx);
midghists = invertselect(midghists, @isempty);
midgf1s = {}; midgoptf1s = [];
for i = 1:length(midghists)
    midgf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), midghists{i});
    midgoptf1s(:,i) = squeeze(max(max(midgf1s{i})));
end
parahists = collect_grating_spikehists(gratingrun, gratingrun.map(get_cell_indices(conerun, {2, 8, 13})), gratingrun.triggers(1:17:end), 'histx', histx);
parahists = invertselect(parahists, @isempty);
paraf1s = {}; paraoptf1s = [];
for i = 1:length(parahists)
    paraf1s{i} = cellfun(@(h)(extract_fourier_components(h, 0.05, 60.35/30, 0.1)), parahists{i});
    paraoptf1s(:,i) = squeeze(max(max(paraf1s{i})));
end


% Plot histogram
figure; hold on;
histedge = 0:0.5:9;
nmidg = histc(midgoptf1s(1,:), histedge);
npara = histc(paraoptf1s(1,:), histedge);
hpara = bar(histedge, npara, 'histc');
hmidg = bar(histedge, nmidg, 'histc');
set(hmidg, 'FaceColor', 'k');
set(hpara, 'FaceColor', 'w');
axis([-0.5 10 0 13]);


%% Raster / PSTH plots (built for 2012-09-13-2)
tpt = 17; % Triggers Per Trial
spatialperiods = unique(gratingrun.stimulus.spatialperiods);

for i = 1:length(gratingrun.offM)
    offM = gratingrun.offM{i};
    if isempty(offM), continue; end
    
    figure
    for j = 1:length(spatialperiods)
        spatialperiod = spatialperiods(j);
        
        for k = 1:4
            spatialphase = spatialperiod * (k-1) / 8;
            
            trials = gratingrun.stimulus.spatialphases == spatialphase & gratingrun.stimulus.spatialperiods == spatialperiod;
            trials = find(trials);
            trials = intersect(trials, gratingrun.stimulus.urgbi{1});
            selected_triggers = repmat((trials-1).*tpt + 1, tpt, 1) + repmat((1:tpt)'-1, 1, length(trials));
            triggers = gratingrun.triggers(selected_triggers(:));
            
            sanesubplot(4, length(spatialperiods), {k j});
            rasterphli(gratingrun, offM, triggers, 'start', -0.1, 'stop', 0.5, 'hist', true, 'hist_bin', 0.05);
            
            if k == 1                
                title(sprintf('spatial period %d', spatialperiod));
            end
            
            if j == 1
                ylabel(sprintf('phase %d/8 period', k-1));
            end
        end
    end
end

for i = 1:length(gratingrun.offP)
    offP = gratingrun.offP{i};
    if isempty(offP), continue; end
    
    figure
    for j = 1:length(spatialperiods)
        spatialperiod = spatialperiods(j);
        
        for k = 1:4
            spatialphase = spatialperiod * (k-1) / 8;
            
            trials = gratingrun.stimulus.spatialphases == spatialphase & gratingrun.stimulus.spatialperiods == spatialperiod;
            trials = find(trials);
            trials = intersect(trials, gratingrun.stimulus.urgbi{1});
            selected_triggers = repmat((trials-1).*tpt + 1, tpt, 1) + repmat((1:tpt)'-1, 1, length(trials));
            triggers = gratingrun.triggers(selected_triggers(:));
            
            sanesubplot(4, length(spatialperiods), {k j});
            rasterphli(gratingrun, offP, triggers, 'start', -0.1, 'stop', 0.5, 'hist', true, 'hist_bin', 0.05);
        end
        
        if k == 1
            title(sprintf('spatial period %d', spatialperiod));
        end
        
        if j == 1
            ylabel(sprintf('phase %d/8 period', k-1));
        end
    end
end


%% F1 components (built for 2012-09-13-2)
tpt = 17; % Triggers Per Trial
spatialperiods = unique(gratingrun.stimulus.spatialperiods);

% Full trials
start = 0;
stop = 8.1;

F = 20;
T = 1/F;
histx = start:T:stop;
L = length(histx);
NFFT = 2^nextpow2(L);
f = F/2 * linspace(0, 1, NFFT/2+1);

screenrefresh = 60.35;
temporalperiod = 30;
f1 = screenrefresh/temporalperiod;
f1x = find(abs(f-f1) < 0.1);

for i = 1:length(gratingrun.offM)
    offM = gratingrun.offM{i};
    if isempty(offM), continue; end

    ax = [];
    figure(offM);
    for j = 1:length(spatialperiods)
        spatialperiod = spatialperiods(j);
        
        for k = 1:4
            spatialphase = spatialperiod * (k-1) / 8;
            
            trials = gratingrun.stimulus.spatialphases == spatialphase & gratingrun.stimulus.spatialperiods == spatialperiod;
            trials = find(trials);
            trials = intersect(trials, gratingrun.stimulus.urgbi{1});
            selected_triggers = (trials-1).*tpt + 1;
            triggers = gratingrun.triggers(selected_triggers(:));
            
            res = rasterphli(gratingrun, offM, triggers, 'start', start, 'stop', stop, 'hist', true, 'hist_bins', histx);
            histy = histc(res(:,1), histx);
            
            histf = fft(histy,NFFT) / L;
            pwr = 2*abs(histf(1:NFFT/2+1));
            baseline = median(pwr);
            midgf1(k,j,i) = sum(pwr(f1x)) - baseline*length(f1x);
            
            ax(k,j) = sanesubplot(4, length(spatialperiods), {k j});
            plot(f(5:end), pwr(5:end));
        end
    end
    
    for j = 1:size(ax,2), title( ax(1,j), sprintf('spatial period %d', spatialperiod)); end
    for k = 1:size(ax,1), ylabel(ax(k,1), sprintf('phase %d/8 period', k-1));           end
    setaxesy(ax, 'tight', true);
end
midgoptf1 = squeeze(max(max(midgf1)));

for i = 1:length(gratingrun.offP)
    offP = gratingrun.offP{i};
    if isempty(offP), continue; end

    ax = [];
    figure(offP);
    for j = 1:length(spatialperiods)
        spatialperiod = spatialperiods(j);
        
        for k = 1:4
            spatialphase = spatialperiod * (k-1) / 8;
            
            trials = gratingrun.stimulus.spatialphases == spatialphase & gratingrun.stimulus.spatialperiods == spatialperiod;
            trials = find(trials);
            trials = intersect(trials, gratingrun.stimulus.urgbi{1});
            selected_triggers = (trials-1).*tpt + 1;
            triggers = gratingrun.triggers(selected_triggers(:));
            
            res = rasterphli(gratingrun, offP, triggers, 'start', start, 'stop', stop, 'hist', true, 'hist_bins', histx);
            histy = histc(res(:,1), histx);
            
            histf = fft(histy,NFFT) / L;
            pwr = 2*abs(histf(1:NFFT/2+1));
            baseline = median(pwr);
            paraf1(k,j,i) = sum(pwr(f1x)) - baseline*length(f1x);
            
            ax(k,j) = sanesubplot(4, length(spatialperiods), {k j});
            plot(f(5:end), pwr(5:end));
        end
    end
    
    for j = 1:size(ax,2), title( ax(1,j), sprintf('spatial period %d', spatialperiod)); end
    for k = 1:size(ax,1), ylabel(ax(k,1), sprintf('phase %d/8 period', k-1));           end
    setaxesy(ax, 'tight', true);
end
paraoptf1 = squeeze(max(max(paraf1)));

figure; hold on;
histedge = 0:0.5:8;
nmidg = histc(midgoptf1(midgoptf1>0), histedge);
npara = histc(paraoptf1(paraoptf1>0), histedge);
hmidg = bar(histedge, nmidg, 'histc');
hpara = bar(histedge, npara, 'histc');
set(hmidg, 'FaceColor', 'k');
set(hpara, 'FaceColor', 'w');
axis([-0.5 8 0 5]);