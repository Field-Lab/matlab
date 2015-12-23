loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'; 'd'};
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};

%% M-
piece = '2011-06-24-6';
datarun = load_data(fullfile(piece, 'data005'), staopts);
rgc = 51;

thresh = 3.5;
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
figure; 
plot_rf(datarun, rgc, 'fit', false, 'autozoom', true, 'scale', 10, 'color_transform', [1 1 1]');

x = 17.5; y = 209.5; del = 32;
axis([x x+del y y+del]);

%% P-
piece = '2012-09-13-2';
datarun = load_data(fullfile(piece, 'data005'), staopts);
rgc = 6031;

thresh = 3.5; % Marks look better at 4, but RF is same
sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
figure; 
plot_rf(datarun, rgc, 'fit', false, 'autozoom', true, 'scale', 10);

x = 181.5; y = 185.5; del = 78;
axis([x x+del y y+del]);

%% M+
piece = '2011-06-24-6';
rgc = 1863;
thresh = 3.5;

d.d03 = load_data([piece '/data003'], staopts);
d.d04 = load_data([piece '/data004'], staopts);
d.d05 = load_data([piece '/data005'], staopts);
d.d03.mapd04 = map_ei(d.d04, d.d03);
d.d05.mapd04 = map_ei(d.d04, d.d05);

d.d04.rgc = 1863;
d.d03.rgc = d.d03.mapd04{get_cell_indices(d.d04, d.d04.rgc)};
d.d05.rgc = d.d05.mapd04{get_cell_indices(d.d04, d.d04.rgc)};

sta = get_sta(d.d03, d.d03.rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d.d03.stas.rfs{d.d03.cell_nums(d.d03.rgc)} = [];
d.d03.stas.marks{d.d03.cell_nums(d.d03.rgc)} = simplemarks;
d.d03.rf = get_rf(d.d03, d.d03.rgc);

sta = get_sta(d.d04, d.d04.rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d.d04.stas.rfs{d.d04.cell_nums(d.d04.rgc)} = [];
d.d04.stas.marks{d.d04.cell_nums(d.d04.rgc)} = simplemarks;
d.d04.rf = get_rf(d.d04, d.d04.rgc);

sta = get_sta(d.d05, d.d05.rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d.d05.stas.rfs{d.d05.cell_nums(d.d05.rgc)} = [];
d.d05.stas.marks{d.d05.cell_nums(d.d05.rgc)} = simplemarks;
d.d05.rf = get_rf(d.d05, d.d05.rgc);

rfcompound = d.d03.rf + mean(d.d04.rf, 3)./3 + mean(d.d05.rf, 3)./3;

datarun = d.d04;
datarun.stas.rfs{get_cell_indices(datarun, d.d04.rgc)} = rfcompound;
% figure;
plot_rf(datarun, d.d04.rgc, 'fit', false, 'autozoom', true, 'scale', 10);

x = 19.5; y = 61.5; del = 39;
axis([x x+del y y+del]);

%% P+
% TODO: Check this one for averaging also
piece = '2011-05-11-6';
d01 = load_data(fullfile(piece, 'data001'), staopts);
d02 = load_data(fullfile(piece, 'data002'), staopts);
d03 = load_data(fullfile(piece, 'data003'), staopts);

d01.mapd02 = map_ei(d02, d01);
d03.mapd02 = map_ei(d02, d03);

d02.onP = 5929;
d01.onP = d01.mapd02{get_cell_indices(d02, d02.onP)};
d03.onP = d03.mapd02{get_cell_indices(d02, d02.onP)};

thresh = 3.5; % Marks look better at 3, but RF looks same
sta = get_sta(d02, d02.onP);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d02.stas.rfs{d02.cell_nums(d02.onP)} = [];
datarun.stas.marks{d02.cell_nums(d02.onP)} = simplemarks;
d02.onPrf = get_rf(d02, d02.onP);

sta = get_sta(d01, d01.onP);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d01.stas.rfs{d01.cell_nums(d01.onP)} = [];
datarun.stas.marks{d01.cell_nums(d01.onP)} = simplemarks;
d01.onPrf = get_rf(d01, d01.onP);

sta = get_sta(d03, d03.onP);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
d03.stas.rfs{d03.cell_nums(d03.onP)} = [];
datarun.stas.marks{d03.cell_nums(d03.onP)} = simplemarks;
d03.onPrf = get_rf(d03, d03.onP);

rfcompound = d01.onPrf + d02.onPrf + d03.onPrf;

datarun = d02;
datarun.stas.rfs{get_cell_indices(datarun, d02.onP)} = rfcompound;
plot_rf(datarun, d02.onP, 'color_transform', [1 1 1]', 'fit', false, 'autozoom', true, 'scale', 10);

x = 55.5; y = 170.5; del = 75;
axis([x x+del y y+del]);

%% Check thresh
thresh = 3.5;

ss = significant_stixels(get_sta(datarun, rgc));
subplot(2,3,1); imagesc(ss); axis equal tight; colormap gray;
autozoom_to_fit(datarun, rgc, 5, 1, 1);

sta = get_sta(datarun, rgc);
%sta = squeeze(sum(get_sta(datarun, rgc), 3)); % Combine channels
maxsta = max(mean(abs(sta), 3), [], 4);
rstd = robust_std(sta(:), 6);
simplemarks = maxsta > (thresh*rstd);

% sta = get_sta(datarun, rgc);
% stacollapse = mean(std(sta, 1, 4), 3);
% rstd = robust_std(stacollapse(:), 6);
% simplemarks = stacollapse > (thresh*rstd);

subplot(2,3,2); imagesc(simplemarks); axis equal tight; colormap gray
autozoom_to_fit(datarun, rgc, 5, 1, 1);
subplot(2,3,3); imagesc(simplemarks); axis equal tight; colormap gray

datarun.stas.marks{datarun.cell_nums(rgc)} = [];
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
rf = get_rf(datarun, rgc);
datarun.stas.rfs{datarun.cell_nums(rgc)} = rf;
subplot(2,3,4); plot_rf(datarun, rgc, 'fit', false, 'autozoom', true);

datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;
datarun.stas.rfs{datarun.cell_nums(rgc)} = [];
simplerf = get_rf(datarun, rgc);
datarun.stas.rfs{datarun.cell_nums(rgc)} = simplerf;
subplot(2,3,5); plot_rf(datarun, rgc, 'fit', false, 'autozoom', true);

if all(rf(:) == simplerf(:)), disp('rf and simplerf are the same!!'); end


%% Old stuff below

%% Searching
%piece = '2011-12-13-2';
%datarunname = '/data004_data007_data008-norefit/data004-from-data004_data007_data008/data004-from-data004_data007_data008';

piece = '2011-05-11-6';
datarunname = '/data002';
datarun = load_data([piece datarunname], staopts);
info(datarun);

%% Plotting
cellspec = {4};
datarun = get_sta_summaries(datarun, cellspec, 'robust_std_method', 5);
plot_rfs(datarun, cellspec, 'overlay', false, 'autozoom', true);


%% Review
piece = '2011-05-11-6';
datarunname = '/data002';
datarun = load_data([piece datarunname], staopts);
for cellid = [1203, 3452, 3841, 5761, 5929, 6361, 6620 394, 3108, 3901, 4127, 5581]
    figure();
    plot_rf(datarun, cellid, 'autozoom', true);
end

%% [1x1!] 2011-08-04-2/data002
% OFF Parasols
% 2041, 2192, 3136
% 1204, 2133

%% 2011-12-13-2/data004cm
% ON Midgets
% 601, 1564, 1683, 2478, 2659, 2778, 2821, 3121, 6541, 6601

% OFF Midgets
% 541, 1006, 1906, 2191, 2656, 2971, 3136, 3586, 3751

% SBC 1715

%% 2012-04-13-1/data002
% ON Parasols
% 4820

% OFF Parasols
% 2176, 4141, 5581, 6093

% OFF Midgets
% 2536
% 2161, 2431, 2778, 3002, 4966, 6046, 6286

% SBC 3184

%% 2012-09-21-2/data007
% OFF Parasol
% 4291, 7141

% ON Midget
% 2147

% OFF Midget
% 541, 1381, 1681, 2221, 2356, 2672, 4096, 4217, 5087, 5491, 7142, 7546

%% 2012-09-13-2/data001
% ON Parasols
% 107, 1322, 2191, 5734, 5748, 6016, 6291, 6482, 6767

% OFF Parasols
% 5807 used for univariance
% 6993 used for rasters
% 451, 1546, 4936, 5041, 5131, 5581, 5761, 5911, 6736

% ON Midgets
% 202, 454, 1743, 2266, 2329, 2492, 2566

% OFF Midgets
% 2117 used for rasters
% 1732, 1741, 2791, 6587, 6796, 7246, 7263, 7459

% SBC
% 6272, 7233

%% 2012-09-13-2/data005
% ON Parasols
% 2177, 5222, 5611, 5733, 6333, 6437, 6856
% 5912

% OFF Parasols
% 5551 Used in univariance fig
% 6992 used in rasters fig
% 5596, 5761, 5926, 6031
% 31, 376, 466, 1411, 1471, 1981, 2118, 4096, 4441, 5041, 5071, 5131, 5371, 7231, 6871
% 46, 286, 2311, 4321, 5401

% ON Midgets
% 1338, 2119, 2388

% OFF Midgets
% 2166 Used for raster
% 47, 151, 452, 1457, 1816, 1966, 6241, 6497, 6796, 6991, 7321, 7576

% SBC
% 6287

% OFF Large
% 7352

%% 2012-09-24-5/data003
% ON Parasols
% 6499, 7237

% OFF Midgets
% 451, 707, 961, 1366, 2596, 2686, 2776, 3151, 4216, 4351, 5056, 7081



%% [RGB] 2011-05-11-6/data002
% ON Parasol
% 1203, 3452, 3841, 5761, 5929, 6361, 6620
% 394, 3108, 3901, 4127, 5581

% OFF Parasol
% 588, 2836, 3451, 3647, 4052, 4172, 4669, 4787, 5672, 6936, 6947

% ON Midget
% 319, 871, 2672, 3048, 3272, 3454, 3707, 3721, 4412, 5147, 5387, 5566, 5671, 5867, 6736, 6754

%% [RGB] 2011-06-24-6/data005
% ON Parasol
% 608, 5898

% OFF Parasol
% 5851

% ON Midget
% 1866

% OFF Midget
% 51, 1293, 5766, 6107, 6438, 6845

%% [RGB] 2011-08-04-2/data005s
% OFF Parasols
% 317

% OFF Midgets
% 452, 586, 1337, 6917, 7188, 7201

% SBC 4939