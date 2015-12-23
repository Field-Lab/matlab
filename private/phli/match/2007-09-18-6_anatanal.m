piece = '2007-09-18-6';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%%
%
d02 = load_data([piece '/data002'], loadopts);
d02a = load_data([piece '/data002-nwpca-all'], loadopts);

%
d02  = calc_csd(d02);
d02a = calc_csd(d02a);

% Load transforms
load(fullfile(server_path(), piece, 'images/TA1'));
% load([server_path '2007-09-18-6/images/TA2.mat'])
load([server_path '2007-09-18-6/images/TA1F1.mat'])
load([server_path '2007-09-18-6/images/TA1F2.mat'])
load([server_path '2007-09-18-6/images/TA1F3.mat'])
load([server_path '2007-09-18-6/images/TA1F4.mat'])
load([server_path '2007-09-18-6/images/TA1F5.mat'])

% Load axons
svg_file_path = fullfile(server_path(), piece, 'images/axon traces.svg/axon traces.svg');
images{1} = {'5749060.tiff',  'f2', maketform('composite',cp2tform(bp_f2,ip_f2,'lwm'), fliptform(TA1)), ip_f2,bp_f2};
images{2} = {'18072020.tiff', 'f3', maketform('composite',cp2tform(bp_f3,ip_f3,'lwm'), fliptform(TA1)), ip_f3,bp_f3};
images{3} = {'5766530.tiff',  'f4', maketform('composite',cp2tform(bp_f4,ip_f4,'lwm'), fliptform(TA1)), ip_f4,bp_f4};
images{4} = {'18000dc0.tiff', 'f5', maketform('composite',cp2tform(bp_f5,ip_f5,'lwm'), fliptform(TA1)), ip_f5,bp_f5};
axons = import_intaglio_axons(svg_file_path, images, TA1);

%% TODO: Update to use d02a after cleaning that up
datarun = d02;
eicorrmatch6 = {};

axid = 173;
cellspecs6{axid} = {3};
matchids6{axid} = 7367;
neighborids6{axid} = [7126 7156 7233 7246 7276 7427 7456 7471 7472 7502 7516 7531 7532 7576 7637];

eicorrmatch6{axid}   = qd_compare_ei_fit(datarun, axons, axid, matchids6{axid}, 'verbose', true);
eicorrnear6{axid}    = qd_compare_ei_fit(datarun, axons, axid, [neighborids6{axid}], 'verbose', true);
eicorrall6{axid}     = qd_compare_ei_fit(datarun, axons, axid, cellspecs6{axid}, 'verbose', true);
eicorrnotnear6{axid} = setdiff(eicorrall6{axid}, [eicorrmatch6{axid}; eicorrnear6{axid}]);

x = 1;
axid = 173;
plot(randn(size(eicorrnotnear6{axid}))/30 + x, eicorrnotnear6{axid}, '.', 'MarkerSize', 3, 'Color', 0.5*[1 1 1]);
plot(x, eicorrnear6{axid}, 'ko', 'MarkerSize', 4);
h = plot(x, eicorrmatch6{axid}, '.', 'MarkerSize', 15, 'Color', [1 0.25 0.25]);


%%
datarun = d02a;
eicorrmatch6a = {};

axid = 173;
cellspecs6a{axid} = {3, 13, 26};
matchids6a{axid} = 7367;
neighborids6a{axid} = [7126 7156 7202 7237 7247 7276 7293 7355 7458 7472 7502 7516 7531 7532 7546 7576];

% find_matches(datarun, axons, axid, 'cell_ids', [matchids6{axid} neighborids6{axid}], 'ei_scale', 0.25)
% cnums = get_cell_indices(datarun, cellspecs6{axid});
% find_matches(datarun, axons, axid, 'cell_ids', datarun.cell_ids(cnums), 'ei_scale', 0.25)

eicorrs = qd_compare_ei_fit(d02a, axons, 173, ax173ids);
[~, corrorder] = sort(eicorrs(:), 1, 'descend');

nplot = length(ax173ids);
nsqrt = sqrt(nplot);
ncol = ceil(nsqrt);
nrow = ceil(nplot / ncol);
figure();
for i = 1:length(ax173ids)
    cid = ax173ids(corrorder(i));
    subplot(nrow, ncol, i);
    hold on
    plot_ei_(get_ei(d02a, cid), d02.ei.position);
    plot(axons{173}(:,1), axons{173}(:,2));
    title(sprintf('%d', cid));
end


%% Predicted ei...
[predax predsom axdists somdists] = generate_predicted_ei(axons{173}, d02a.ei.position);
plot_ei_(predax+predsom, d02.ei.position);
plot_ei_(predax, d02.ei.position);
plot_ei_(predsom, d02.ei.position);