piece = '2007-09-18-4';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%%
d02a = load_data([piece '/data002-nwpca-duplicates/data002/data002'], staopts);
d02a = restore_stacks(d02a);
[~, ai] = load_array_info(d02a);
% d02a = load_txt_cell_types(d02a, [d02a.names.rrs_prefix '.classification.PHL.txt']);

% Load transforms
load(fullfile(server_path(), piece, 'images/TA1.mat'));
load(fullfile(server_path(), piece, 'images/TA2.mat'));
load(fullfile(server_path(), piece, 'images/TA1F1.mat'));
load(fullfile(server_path(), piece, 'images/TA1F2.mat'));
load(fullfile(server_path(), piece, 'images/TA1F3.mat'));
load(fullfile(server_path(), piece, 'images/TA1F4.mat'));
load(fullfile(server_path(), piece, 'images/TA1F5.mat'));
load(fullfile(server_path(), piece, 'images/TA1F6.mat'));

% Load axons
svg_file_path = fullfile(server_path(), piece, 'images/axon traces.svg/axon traces.svg');
images{1} = {'1ff3aff0.tiff','f2',maketform('composite',cp2tform(bp_f2,ip_f2,'lwm'),fliptform(TA1)),ip_f2,bp_f2};
images{2} = {'2345ade0.tiff','f3',maketform('composite',cp2tform(bp_f3,ip_f3,'lwm'),fliptform(TA1)),ip_f3,bp_f3};
images{3} = {'1b9c2740.tiff','f4',maketform('composite',cp2tform(bp_f4,ip_f4,'lwm'),fliptform(TA1)),ip_f4,bp_f4};
images{4} = {'1588a450.tiff','f5',maketform('composite',cp2tform(bp_f5,ip_f5,'lwm'),fliptform(TA1)),ip_f5,bp_f5};
images{5} = {'20bcdd10.tiff','f6',maketform('composite',cp2tform(bp_f6,ip_f6,'lwm'),fliptform(TA1)),ip_f6,bp_f6};
axons = import_intaglio_axons(svg_file_path,images,TA1);
% fix a booboo
axons{14} = axons{14}([1 3:end],:);


%% Depth plot
axid = 9;
somdenfn = ['somden' num2str(axid)];
somdenpath = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '.tif']);
somdenpath_despeck = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '_despeck.tif']);

if exist(somdenpath_despeck, 'file')
    somdenpath = somdenpath_despeck;
end

clear somden
for i = 1:100
    somden(:,:,i) = imread(somdenpath, i);
end
somden = reshape(somden, [447 408 2 50]);
somden(:,:,3,:) = somden(:,:,2,:);
somden(:,:,2,:) = somden(:,:,1,:);
somden(:,:,1,:) = 0;

clear zgval zbval
for i = 1:50
    zgval(i) = sum(sum(somden(:,:,2,i)));
    zbval(i) = sum(sum(somden(:,:,3,i)));
end
zgval = zgval - min(zgval);
zgval = zgval ./ max(zgval);
zbval = zbval - min(zbval);
zbval = zbval ./ max(zbval);

subplot(2,1,1);
plot(1:50, zbval);
hold on
plot(1:50, zgval, 'g');
title('ON Midget');


%% Depth plot 2
axid = 17;
somdenfn = ['somden' num2str(axid)];
somdenpath = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '.tif']);
somdenpath_despeck = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '_despeck.tif']);

if exist(somdenpath_despeck, 'file')
    somdenpath = somdenpath_despeck;
end

clear somden
for i = 1:88
    somden(:,:,i) = imread(somdenpath, i);
end
somden = reshape(somden, size(somden,1), size(somden,2), 2, []);
somden(:,:,3,:) = somden(:,:,2,:);
somden(:,:,2,:) = somden(:,:,1,:);
somden(:,:,1,:) = 0;

clear zgval zbval
for i = 1:size(somden,4)
    zgval(i) = sum(sum(somden(:,:,2,i)));
    zbval(i) = sum(sum(somden(:,:,3,i)));
end
zgval = zgval - min(zgval);
zgval = zgval ./ max(zgval);
zbval = zbval - min(zbval);
zbval = zbval ./ max(zbval);

subplot(2,1,2);
plot(1:length(zbval), zbval);
hold on
plot(1:length(zgval), zgval, 'g');
title('OFF Midget');


%% ON Midget figure
datarun = d02a;
cellspec = {3};
axid = 9;
cid = 2296;
neighborids{axid} = [2314 2312 2197 2311 2431];
neighborids2{axid} = [2181 2191 2086 2236 2467 2551 2536 2417 2401 2282];

datarun = load_ei(datarun, cellspec);
datarun = calc_csd(datarun, 'cellspec', cellspec);

% find_matches(datarun, axons, axid, 'cell_ids', [cid neighborids{axid} neighborids2{axid}], 'ei_scale', 0.25)
eicorrmatch{axid} = qd_compare_ei_fit(datarun, axons, axid, cid, 'verbose', true);
eicorrnear{axid} = qd_compare_ei_fit(datarun, axons, axid, [neighborids{axid} neighborids2{axid}], 'verbose', true);
eicorrall{axid}  = qd_compare_ei_fit(datarun, axons, axid, cellspec, 'verbose', true);
eicorrnotnear{axid} = setdiff(eicorrall{axid}, [eicorrmatch{axid}; eicorrnear{axid}]);


% Plot
f = figure();
set(f, 'Color', 'w');

subplot(1,3,2);
plot_rf_summaries(datarun, cid, 'label', true, 'fill_color', [1 0.75 0.75]); 
plot_rf_summaries(datarun, [neighborids{axid} neighborids2{axid}], 'label', true); 
axis([54.3322   63.3843    7.1534   16.2055]);
axis off;
view(160,-90)
% plot_rf_summaries(datarun, cellspec, 'label', true); 
% autozoom_to_fit(datarun, cid, 15);

% Matches
eibox = [-1100 600 -600 100];
eiview = {180, -90};
sanesubplot(3,3,{1 1});
cla
plot_ei(datarun, cid, 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', cid, eicorrmatch{axid}));

idx = 4;
sanesubplot(3,3,{2 1});
cla
plot_ei(datarun, neighborids2{axid}(idx), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(idx), eicorrnear{axid}(5+idx)));

idx = 3;
sanesubplot(3,3,{3 1});
cla
plot_ei(datarun, neighborids2{axid}(idx), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(idx), eicorrnear{axid}(5+idx)));

idx = 7;
sanesubplot(3,3,{1 3});
cla
plot_ei(datarun, neighborids2{axid}(idx), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(idx), eicorrnear{axid}(5+idx)));

idx = 8;
sanesubplot(3,3,{2 3});
cla
plot_ei(datarun, neighborids2{axid}(idx), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(idx), eicorrnear{axid}(5+idx)));

idx = 1;
sanesubplot(3,3,{3 3});
cla
plot_ei(datarun, neighborids2{axid}(idx), 'scale', 0.75, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
axis(eibox);
view(eiview{:});
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(idx), eicorrnear{axid}(5+idx)));


%% ON Midget figure 2
datarun = d02a;
cellspec = {3};
axid = 11;
cid = 2913;
neighborids{axid} = [2912 3032 3016 2897];
neighborids2{axid} = [2926 2942 3170 3271 3151 2882 2776 2794];

datarun = load_ei(datarun, cellspec);
datarun = calc_csd(datarun, 'cellspec', cellspec);

% find_matches(datarun, axons, axid, 'cell_ids', [cid neighborids{axid} neighborids2{axid}], 'ei_scale', 0.25)
% find_matches(datarun, axons, axid, 'cell_ids', datarun.cell_types{cellspec{:}}.cell_ids, 'ei_scale', 0.25)
eicorrmatch{axid} = qd_compare_ei_fit(datarun, axons, axid, cid, 'verbose', true);
eicorrnear{axid} = qd_compare_ei_fit(datarun, axons, axid, [neighborids{axid} neighborids2{axid}], 'verbose', true);
eicorrall{axid}  = qd_compare_ei_fit(datarun, axons, axid, cellspec, 'verbose', true);
eicorrnotnear{axid} = setdiff(eicorrall{axid}, [eicorrmatch{axid}; eicorrnear{axid}]);


% Plot
f = figure();
set(f, 'Color', 'w');

% RF Mosaic
sanesubplot(1,3,{1 2});
plot_rf_summaries(datarun, cid, 'label', true, 'fill_color', [1 0.75 0.75]); 
plot_rf_summaries(datarun, [neighborids{axid} neighborids2{axid}], 'label', true); 
axis([54.4299   61.8949   17.6786   25.1436]);
axis off;
view(180,90)
% plot_rf_summaries(datarun, cellspec, 'label', true); 
% autozoom_to_fit(datarun, cid, 15);

% Matches
sanesubplot(3,3,{1 1});
plot_ei(datarun, cid, 'scale', 1.25, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', cid, eicorrmatch{axid}));

sanesubplot(3,3,{2 1});
ind = 6;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(4+ind)));

sanesubplot(3,3,{3 1});
ind = 7;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(4+ind)));

sanesubplot(3,3,{1 3});
cla
plot_ei(datarun, neighborids2{axid}(3), 'scale', 1.25, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(3), eicorrnear{axid}(4+3)));

sanesubplot(3,3,{2 3});
cla
plot_ei(datarun, neighborids{axid}(2), 'scale', 1.25, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids{axid}(2), eicorrnear{axid}(2)));

sanesubplot(3,3,{3 3});
ind = 1;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(4+ind)));


%% OFF Midget figure
datarun = d02a;
cellspec = {4};
axid = 17;
cid = 3559;
% neighborids{axid} = [3301, 3409, 3413, 3456, 3531, 3545, 3556, 3576, 3681];
neighborids{axid} = [3456 3556 3531 3545];
neighborids2{axid} = [3301 3413 3409 3681 3576 3338 3288];
% neighborids2{axid} = [3275, 3288, 3338, 3412, 3588, 3590, 3635, 3682];

datarun = load_ei(datarun, {4});
datarun = calc_csd(datarun, 'cellspec', {4});

% find_matches(datarun, axons, axid, 'cell_ids', [cid neighborids{axid} neighborids2{axid}], 'ei_scale', 0.25)
eicorrmatch{axid} = qd_compare_ei_fit(datarun, axons, axid, cid, 'verbose', true); 
eicorrnear{axid} = qd_compare_ei_fit(datarun, axons, axid, [neighborids{axid} neighborids2{axid}], 'verbose', true);
eicorrall{axid}  = qd_compare_ei_fit(datarun, axons, axid, {4}, 'verbose', true);
eicorrnotnear{axid} = setdiff(eicorrall{axid}, [eicorrmatch{axid}; eicorrnear{axid}]);


% Plot
f = figure();
set(f, 'Color', 'w');

% RF Mosaic
sanesubplot(1,3,{1 2});
plot_rf_summaries(datarun, cid, 'label', true, 'fill_color', [1 0.75 0.75]); 
plot_rf_summaries(datarun, [neighborids{axid} neighborids2{axid}], 'label', true); 
axis([51.0687   57.8265   23.9048   30.9417]);
axis off;
view(150,90);
% plot_rf_summaries(datarun, cellspec, 'label', true); 
% autozoom_to_fit(datarun, cid, 15);

% Matches
sanesubplot(3,3,{1 1});
cla
plot_ei(datarun, cid, 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', cid, eicorrmatch{axid}));

sanesubplot(3,3,{1 3});
cla
ind = 2;
plot_ei(datarun, neighborids{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids{axid}(ind), eicorrnear{axid}(ind)));

sanesubplot(3,3,{2 3});
cla
ind = 1;
plot_ei(datarun, neighborids{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids{axid}(ind), eicorrnear{axid}(ind)));

sanesubplot(3,3,{3 3});
cla
ind = 1;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(length(neighborids{axid}) + ind)));

sanesubplot(3,3,{2 1});
cla
ind = 3;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(length(neighborids{axid}) + ind)));

sanesubplot(3,3,{3 1});
cla
ind = 7;
plot_ei(datarun, neighborids2{axid}(ind), 'scale', 1, 'coordinates', 'array', ...
    'neg_color', 0.5 * [1 1 1], 'pos_color', 0.5 * [1 1 1], ...
    'axon', axons{axid}, 'axon_color', 'r');
axis off
title(sprintf('ON Midget %d\nCorrelation: %.2f', neighborids2{axid}(ind), eicorrnear{axid}(length(neighborids{axid}) + ind)));


%% 1D population plot
datarun = d02a;
corrfunc = [];
eicorrmatch = {};

% Search for more?
% figure();
% plot_rf_summaries(datarun, cellspecs{axid}, 'label', true); 
% plot_rf_summaries(datarun, matchids{axid}, 'fit_color', 'b', 'fit_width', 2); 
% autozoom_to_fit(datarun, matchids{axid}, 15);
% 
% find_matches(datarun, axons, axid, 'cell_ids', [matchids{axid} neighborids{axid} neighborids2{axid}], 'ei_scale', 0.25)
% find_matches(datarun, axons, axid, 'cell_ids', datarun.cell_types{cellspecs{axid}{1}}.cell_ids, 'ei_scale', 0.25)


% Copied from above
axid = 9;
cellspecs{axid} = {3};
matchids{axid} = 2296;
neighborids{axid} = [2314 2312 2197 2311 2431];
neighborids2{axid} = [2181 2191 2086 2236 2467 2551 2536 2417 2401 2282];

axid = 11;
cellspecs{axid} = {3};
matchids{axid} = 2913;
neighborids{axid} = [2912 3032 3016 2897];
neighborids2{axid} = [2926 2942 3170 3271 3151 2882 2776 2794];

axid = 17;
cellspecs{axid} = {4};
matchids{axid} = 3559;
neighborids{axid} = [3456 3556 3531 3545];
neighborids2{axid} = [3301 3413 3409 3681 3576 3338 3288];

% Additional 1D population data

% TODO: Look into correlations
axid = 15;
cellspecs{axid} = {3};
matchids{axid} = 2342;
neighborids{axid} = [2492 2602 2581 2686 2566 2467];
neighborids2{axid} = [2357 2491 2614 2807 2671 2673 2552 2551 2311 2236];

% TODO: Looks decent, but correlations aren't good; look into
axid = 80;
cellspecs{axid} = {4};
matchids{axid} = 413;
neighborids{axid} = [411 574 271 214 289 348];
neighborids2{axid} = [171 649 692 158 196];

axid = 44;
cellspecs{axid} = {4};
matchids{axid} = 4446;
neighborids{axid} = [4400 4326 4402];
neighborids2{axid} = [4517 4336 4281 4208 4638 4640];

axid = 14;
cellspecs{axid} = {14};
matchids{axid} = 2225;
neighborids{axid} = [2026 2388 2584];
neighborids2{axid} = [1773 1473 2871 2806 2881];

% TODO: Looks good, but correlation needs to be looked into
axid = 41;
cellspecs{axid} = {14};
matchids{axid} = 1326;
neighborids{axid} = [4351 4967 785 977 1219];
neighborids2{axid} = [4397 4816 4279 2871 2388 1473 677 483 5271 5573];

for axid = 1:length(matchids)
    if length(eicorrmatch) >= axid && ~isempty(eicorrmatch{axid}), continue; end
    
    m = matchids{axid};
    if isempty(m), continue; end
    axid
    
    n = neighborids{axid};
    n2 = neighborids2{axid};
    cs = cellspecs{axid};
    eicorrmatch{axid} = qd_compare_ei_fit(datarun, axons, axid, m, 'matchopts', {'corrfunc', corrfunc});
    eicorrnear{axid} = qd_compare_ei_fit(datarun, axons, axid, [neighborids{axid} neighborids2{axid}], 'matchopts', {'corrfunc', corrfunc});
    eicorrall{axid}  = qd_compare_ei_fit(datarun, axons, axid, cs, 'matchopts', {'corrfunc', corrfunc});
    eicorrnotnear{axid} = setdiff(eicorrall{axid}, [eicorrmatch{axid}; eicorrnear{axid}]);
end


% Plot 1D
% (Leaving spaces for parasols from 207-09-18-0)
axonids = [9    11    15    17   44   80    14    41  ];
xs =      [0.5   1.5   2     3   3.5   4     7     7.5];

f1d = figure();
axis([0 8 0 1]);
set(gca(), 'XTick', []);
set(gca(), 'Color', 'none');
set(gca(), 'YTick', -0.2:0.2:0.9);
p = get(f1d, 'Position');
set(f1d, 'Position', [p(1:2) 750 350])
set(f1d, 'Color', 'w');
hold on

for i = 1:length(axonids)
    axid = axonids(i);
    x = xs(i);
    plot(randn(size(eicorrnotnear{axid}))/30 + x, eicorrnotnear{axid}, '.', 'MarkerSize', 3, 'Color', 0.5*[1 1 1]);
    plot(x, eicorrnear{axid}, 'ko', 'MarkerSize', 4);
    h = plot(x, eicorrmatch{axid}, '.', 'MarkerSize', 15, 'Color', [1 0.25 0.25]);
end

textheight = 0.98;
text(1.25, textheight, 'ON Midget', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(xs([1 3]), (textheight - 0.02)*[1 1], 'k--', 'Clipping', 'off');
text(xs(5), textheight, 'OFF Midget', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(xs([4 6]), (textheight - 0.02)*[1 1], 'k--', 'Clipping', 'off');
text(mean(xs([7 8])), textheight, 'SBC', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(xs([7 8]), (textheight - 0.02)*[1 1], 'k--', 'Clipping', 'off');


%% Histogram
% histx = -0.225:0.025:0.8;
% nall = hist(eicorrall{axid}, histx);
% nnear = hist(eicorrnear{axid}, histx);
% nnotnear = nall - nnear;
% bar(histx, nall, 'hist')
% hold on;
% b = bar(histx, nnear, 'hist');
% set(b, 'FaceColor', 'w');
% axis tight
% set(gca, 'YLim', [0 60]);


%% Merge classifications
% Done

% offM = cell2mat(load_txt_cell_types([d02a.names.rrs_prefix '.classification.offM.PHL.txt'], 'cell_type_depth', 5, 'order_cell_types', false));
% onM  = cell2mat(load_txt_cell_types([d02a.names.rrs_prefix '.classification.onM.PHL.txt'],  'cell_type_depth', 5, 'order_cell_types', false));
% 
% allids = [offM.cell_ids];
% offMids = offM(cellfun(@(n) isequal(n, 'OFF midget clean '), {offM.name})).cell_ids;
% onMids  = onM( cellfun(@(n) isequal(n, 'ON midget clean '),  {onM.name})).cell_ids;
% unclassids = setdiff(allids, [offMids onMids]);
% 
% offMclass = struct('name', 'midget', 'cells', offMids, 'subclasses', struct([]));
% offclass = struct('name', 'OFF', 'cells', [], 'subclasses', offMclass);
% onMclass = struct('name', 'midget', 'cells', onMids, 'subclasses', struct([]));
% onclass = struct('name', 'ON', 'cells', [], 'subclasses', onMclass);
% unclass = struct('name', 'unclassified', 'cells', unclassids, 'subclasses', struct([]));
% class = struct('name', 'All', 'cells', [], 'subclasses', [onclass offclass unclass]);
% 
% write_txt_classification(class, [d02a.names.rrs_prefix '.classification.PHL.txt'], 'verbose', true)


%% Find electrodes to cell-finder, hopefully find more cells
cids = [2296 3559];
elecs = arrayfun(...
    @(cid) edu.ucsc.neurobiology.vision.io.NeuronFile.getNeuronIDElectrode(cid), ...
    cids);

radius = 2;
adjelecs = arrayfun(...
    @(e) ai.electrodeMap.getAdjacentsTo(e, radius), ...
    elecs, 'UniformOutput', false);

adjelecs = unique(cell2mat(adjelecs(:)));
