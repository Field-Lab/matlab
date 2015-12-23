%% Basics
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};


%% Most double cone have only max contrast :(
% Includes:
%   2012-09-13-2
%   2012-09-24-5


%% Some hope in the early data where we were just doing additivity:
%   2011-05-11-6
%       d02/d07     ON P 18 19 20 22 have /some/ response...
%
%   2011-07-05-4
%       d01,d03/d05     Nothing much
%
%   2011-07-14-6
%       d00,d02/d04     Nothing much
%
%   2011-08-04-2
%       Worth checking, but stability not great
%
%   2011-08-04-6
%       Worth checking, but responses not great
%

%% 2012-09-24-5
piece = '2012-09-24-5';
d.d00s = load_data(fullfile(piece, 'streamed/data000/data000'), staopts);
d.d01 = load_data(fullfile(piece, 'data001'), staopts);
d.d03 = load_data(fullfile(piece, 'data003'), staopts);
d.d05 = load_data(fullfile(piece, 'data005'), staopts);
d.d03_06_07 = load_data(fullfile(piece, 'd03-06-07-norefit/d03-06-07-norefit'), loadopts);

d.d03_06_07.mapd03 = map_ei(d.d03, d.d03_06_07);
d.d03_06_07.mapd05 = map_ei(d.d05, d.d03_06_07);

d.d06cm = load_data(fullfile(piece, 'd03-06-07-norefit/data006/data006'), loadopts);
d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);

d.d03.onMs = [367 1011 2614 2855 3064 4597];
d.d03.onPs = [3739 7237];

d.d06cm.onMs = [d.d03_06_07.mapd03{get_cell_indices(d.d03, d.d03.onMs)}];
d.d06cm.onPs = [d.d03_06_07.mapd03{get_cell_indices(d.d03, d.d03.onPs)}];

urgb = d.d06cm.stimulus.urgb;
doubleups = find(urgb.doubles & any(urgb.incr) & ~urgb.uds);
titles = {'3+4' '2+4' '2+3' '1+4' '1+3' '1+2'};
for i = 1:length(d.d06cm.onMs)
    onM = d.d06cm.onMs(i);
    figure(onM);
    urgb_raster_subplot(d.d06cm, onM, d.d06cm.triggers(1:2:end), num2cell(doubleups), 'titles', titles);
end
for i = 1:length(d.d06cm.onPs)
    onP = d.d06cm.onPs(i);
    figure(onP);
    urgb_raster_subplot(d.d06cm, onP, d.d06cm.triggers(1:2:end), num2cell(doubleups), 'titles', titles);
end
% Essentially, the single good onP comes in at ~12 Hz peak, which is right
% at the bottom of the distribution of onM, which range from around 10-35
% Hz.


d.d05.onMs = [228 2462 2853 3064 3738 4531 6346];


d.d01.onPs = [336]; % For d.d04 rasters
d.d01.onMs = [2977 3933]; % For d.d04 rasters; should look for more


%% 2011-07-14-6
% piece = '2011-07-14-6';
% d.d00 = load_data(fullfile(piece, 'data000'), staopts);
% d.d02 = load_data(fullfile(piece, 'streamed/data002/data002'), staopts);
% d.d04 = load_data(fullfile(piece, 'data004'), loadopts);
% d.d04 = read_stim_lisp_output(d.d04);
% d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
% 
% d.d04.mapd00 = map_ei(d.d00, d.d04);
% additivity_raster_plot(d.d04, d.d04.stimulus, [d.d04.mapd00(get_cell_indices(d.d00, {1}))], 'triggers', d.d04.triggers(1:3:end), 'figure_by_cell_id', true);
% % One ON Parasol with surround responses...
% % One ON Midget with surround responses...
% 
% d.d04.mapd02 = map_ei(d.d02, d.d04);
% additivity_raster_plot(d.d04, d.d04.stimulus, [d.d04.mapd02(get_cell_indices(d.d02, {3}))], 'triggers', d.d04.triggers(1:3:end), 'figure_by_cell_id', true);
% % Nothing much for ONs

%% 2011-07-05-4
% piece = '2011-07-05-4';
% d.d01 = load_data(fullfile(piece, 'streamed/data001/data001'), staopts);
% d.d03 = load_data(fullfile(piece, 'data003'), staopts);
% d.d05 = load_data(fullfile(piece, 'data005'), loadopts);
% 
% d.d05 = read_stim_lisp_output(d.d05);
% d.d05.stimulus = parse_stim_rgbs(d.d05.stimulus);
% 
% % Searching
% d.d05.mapd03 = map_ei(d.d03, d.d05);
% additivity_raster_plot(d.d05, d.d05.stimulus, [d.d05.mapd03(get_cell_indices(d.d03, {22}))], 'figure_by_cell_id', true);
% % Not much for ON Parasol, maybe a few surround responses
% % Some center and surrounds for OFF Parasol
% % Nothing for ON Midget
% % Nice centers and surrounds for OFF Midget


%% 2011-06-24-6 
piece = '2011-06-24-6';
d.d06 = load_data(fullfile(piece, 'data006'), loadopts);
d.d06 = read_stim_lisp_output(d.d06);
d.d06.stimulus = parse_stim_rgbs(d.d06.stimulus);

% Searching
%d.d06.mapd02 = map_ei(d02, d.d06);
%additivity_raster_plot(d.d06, d.d06.stimulus, [d06.mapd02{get_cell_indices(d02, {1})}], 'figure_by_cell_id', true);

%% 2011-05-11-6 data007
piece = '2011-05-11-6';
d.d07 = load_data(fullfile(piece, 'data007'), loadopts);
d.d07 = read_stim_lisp_output(d.d07, 'add_1_2_data002_600x600');

% Pick timescale several factors smaller than mode ISI
delT = 0.005;

% Boxcar template
% We're going to use Gaussian smoothed peaks for this version, so pick fine time binning.
boxstart = 0.05;
boxend = 0.25;
template.x = -0.1:delT:0.5;
template.y = zeros(size(template.x));
template.y(template.x > boxstart & template.x < boxend) = 1;
template.y = template.y ./ sum(template.y);

crsx = [0.5 0.25 0.125 0 -0.125 -0.25 -0.5];
onPs = [4381 4489 5071 5776];
for i = 1:length(onPs)
    onP = onPs(i);
    [~, ~, histy, res] = additivity_raster_plot(d.d07, d.d07.stimulus, onP, 'hist_bins', template.x, 'plot', false);
    tempcrs = NaN(size(histy));
    for j = 1:numel(histy)
        tempcrs(j) = template.y(:)' * histy{j}(:);
    end
    paracrs{i} = tempcrs([15 16 17 4 19 20 21]); % Kludge to get doubles plus blank
    
    % For Gaussian filter peak calc
    parahists{i} = histy;
    
%    parafits(i,:) = polyfit(crsx(1:4), paracrs{i}(1:4) - paracrs{i}(4), 1);
end

onMs = [4652 5476 5866 6421 6511 4818 6721 3272 7458 7486 1397];
for i = 1:length(onMs)
    onM = onMs(i);
    [~, ~, histy, res] = additivity_raster_plot(d.d07, d.d07.stimulus, onM, 'hist_bins', template.x, 'plot', false);
    tempcrs = NaN(size(histy));
    for j = 1:numel(histy)
        tempcrs(j) = template.y(:)' * histy{j}(:);
    end
    midgcrs{i} = tempcrs([15 16 17 4 19 20 21]); % Kludge to get doubles plus blank
 
    % For Gaussian filter peak calc
    midghists{i} = histy;
    
%    midgfits(i,:) = polyfit(crsx(1:4), midgcrs{i}(1:4) - midgcrs{i}(4), 1);
end


% Gaussian filter and get peaks
gaussx = -0.5:delT:0.5;
gaussy = normpdf(gaussx, 0, 0.03); % sigma close to ISI calculated above
filtx = convx(template.x, gaussx);
for i = 1:length(onPs)
    filtparahists{i} = cellfun(@(H)(delT .* conv(H, gaussy)), parahists{i}, 'UniformOutput', false);
end
for i = 1:length(onMs)
    filtmidghists{i} = cellfun(@(H)(delT .* conv(H, gaussy)), midghists{i}, 'UniformOutput', false);
end
stimx = false(size(filtx));
stimx(filtx > 0 & filtx < 0.25) = true;
for i = 1:length(onPs)
    parapeaks{i} = cellfun(@(H)(max(H(stimx))), filtparahists{i});
end
for i = 1:length(onMs)
    midgpeaks{i} = cellfun(@(H)(max(H(stimx))), filtmidghists{i});
end

% Kludge to get just positive contrast doubles
parapeaks = cell2mat(cellfun(@(P)(P([15 16 17])), parapeaks', 'UniformOutput', false));
midgpeaks = cell2mat(cellfun(@(P)(P([15 16 17])), midgpeaks', 'UniformOutput', false));


figure;
hold on;
for i = 1:3 % Cut last parasol
    plot(crsx(1:3), parapeaks(i,:) - mean(parahists{i}{4}), 'k^--');
end
for i = 1:length(midgcrs)
    plot(crsx(1:3), midgpeaks(i,:) - mean(midghists{i}{4}), 'ko-');
end


%% Old boxcar (not peaks) C/R of ON Parasol doubles versus ON Midget doubles (2011-05-11-6)
boxstart = 0.05;
boxend = 0.25;
template.x = -0.1:0.05:0.5;
template.y = zeros(size(template.x));
template.y(template.x > boxstart & template.x < boxend) = 1;
template.y = template.y ./ sum(template.y);

crsx = [0.5 0.25 0.125 0 -0.125 -0.25 -0.5];
onPs = [4381 4489 5071 5776];
for i = 1:length(onPs)
    onP = onPs(i);
    [~, ~, histy, res] = additivity_raster_plot(d07, s07, onP, 'hist_bins', template.x, 'plot', false);
    tempcrs = NaN(size(histy));
    for j = 1:numel(histy)
        tempcrs(j) = template.y(:)' * histy{j}(:);
    end
    paracrs{i} = tempcrs([15 16 17 4 19 20 21]); % Kludge to get doubles plus blank
    
    parafits(i,:) = polyfit(crsx(1:4), paracrs{i}(1:4) - paracrs{i}(4), 1);
end

onMs = [4652 5476 5866 6421 6511 4818 6721 3272 7458 7486 1397];
for i = 1:length(onMs)
    onM = onMs(i);
    [~, ~, histy, res] = additivity_raster_plot(d07, s07, onM, 'hist_bins', template.x, 'plot', false);
    tempcrs = NaN(size(histy));
    for j = 1:numel(histy)
        tempcrs(j) = template.y(:)' * histy{j}(:);
    end
    midgcrs{i} = tempcrs([15 16 17 4 19 20 21]); % Kludge to get doubles plus blank
    
    midgfits(i,:) = polyfit(crsx(1:4), midgcrs{i}(1:4) - midgcrs{i}(4), 1);
end

figure;
hold on;
for i = 1:length(paracrs)
    plot(crsx(1:3), paracrs{i}(1:3) - paracrs{i}(4), 'k^--');
end
for i = 1:length(midgcrs)
    plot(crsx(1:3), midgcrs{i}(1:3) - midgcrs{i}(4), 'ko-');
end