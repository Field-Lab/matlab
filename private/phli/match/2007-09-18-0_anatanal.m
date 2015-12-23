piece = '2007-09-18-0';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%%
d01a = load_data([piece '/data001-nwpca-all'], staopts);
[~, ai] = load_array_info(d01a);

% Transforms
load([server_path '2007-09-18-0/images/TA1.mat'])
%load([server_path '2007-09-18-0/images/TA2.mat'])
%load([server_path '2007-09-18-0/images/TA1F1.mat'])
load([server_path '2007-09-18-0/images/TA1F2.mat'])
load([server_path '2007-09-18-0/images/TA1F3.mat'])
load([server_path '2007-09-18-0/images/TA1F4.mat'])
load([server_path '2007-09-18-0/images/TA1F5.mat'])
load([server_path '2007-09-18-0/images/TA1F6.mat'])

% Axons
svg_file_path = fullfile(server_path(), piece, 'images/axon traces.svg/axon traces.svg');
images{1} = {'18f20030.tiff','f2',maketform('composite',cp2tform(bp_f2,ip_f2,'lwm',9),fliptform(TA1)),ip_f2,bp_f2};
images{2} = {'16726d00.tiff','f3',maketform('composite',cp2tform(bp_f3,ip_f3,'lwm'),fliptform(TA1)),ip_f3,bp_f3};
images{3} = {'21331ef0.tiff','f4',maketform('composite',cp2tform(bp_f4,ip_f4,'lwm'),fliptform(TA1)),ip_f4,bp_f4};
images{4} = {'167eeaa0.tiff','f5',maketform('composite',cp2tform(bp_f5,ip_f5,'lwm'),fliptform(TA1)),ip_f5,bp_f5};
images{5} = {'18f20ab0.tiff','f6',maketform('composite',cp2tform(bp_f6,ip_f6,'lwm'),fliptform(TA1)),ip_f6,bp_f6};
[axons axT] = import_intaglio_axons(svg_file_path,images,TA1);

% I extended this axon, but did it in Inkscape instead of Intaglio so becomes pain in ass to bring in changes
% Relative coordinates from InkScape
ax164 = cumsum([270.507,388.32;1.203,-1.916;1.264,-1.815;1.174,-1.406;0.579,-0.591;0.979,-0.612;0.204,-1.325;0.835,-1.774;1.061,-2.018;0.998,-1.345;1.489,-1.489;1.467,-1.345;2.263,-1.753;1.754,-1.509;2.079,-1.651;1.672,-2.263;2.018,-2.406;1.121,-1.57;1.142,-1.304;1.182,-1.468;1.325,-2.671;1.325,-2.161;0.653,-0.958;1.61,-1.529;1.162,-1.386;1.815,-2.263;0.917,-1.488;0.653,-1.57;0.611,-1.08;0.673,-1.529;0.448,-1.285;0.225,-1.569;0.55,-1.163;1.182,-2.364;0.755,-1.468;1.039,-1.509;0.816,-0.795;2.079,-2.222;0.959,-0.917;1.427,-1.57;0.856,-0.877;0.815,-1.488;0.673,-1.223;0.347,-1.346;0.897,-1.488;0.733,-1.264;1.32298,-2.01865;0.60234,-1.0768;0.31214,-0.4813;0.82598,-1.64376;0.73252,-2.06196;0.72546,-1.57271;1.23993,-2.36689;0.14665,-1.01319;0.34373,-0.63681;0.33421,-0.94561;0.38258,-0.99979;0.0317,-0.91356;0.15413,-0.86708;0.24418,-1.36492;0.16874,-0.89474;0.19231,-0.99212;0.0874,-0.61356;0.19009,-0.94295;0.22805,-0.99113;0.40371,-1.00434;0.406,-0.86377;0.5128,-1.30727;0.43033,-1.22201;0.33781,-1.20446;0.11599,-1.45477;0.45187,-1.57559;0.24627,-1.65681;0.14549,-1.11491;0.16518,-1.09745;0.2013,-0.96534;0.18143,-0.75449;0.18504,-1.80019;0.14456,-1.44272;0.18667,-1.13963;0.16918,-0.9402;0.24538,-1.37531;0.33843,-1.14492;0.31441,-1.16556;0.34997,-1.38305;0.36163,-1.38061;0.36793,-1.98965;0.43527,-1.50767;0.40491,-1.44456;0.55876,-1.57002;0.66855,-1.81747;1.16539,-3.08136]);
ax164T = tforminv(axT, ax164);
axons{164} = ax164T;

%% 1D population plot
datarun = d01a;
eicorrmatch0 = {};

% Look for more?
% figure();
% plot_rf_summaries(datarun, cellspecs0{axid}, 'label', true); 
% plot_rf_summaries(datarun, matchids0{axid}, 'fit_color', 'b', 'fit_width', 2); 
% autozoom_to_fit(datarun, matchids0{axid}, 15);
% 
% find_matches(datarun, axons, axid, 'cell_ids', [matchids0{axid} neighborids0{axid} neighborids2_0{axid}], 'ei_scale', 0.25)
% find_matches(datarun, axons, axid, 'cell_ids', datarun.cell_types{cellspecs0{axid}{1}}.cell_ids, 'ei_scale', 0.25)


% TODO: Look into; not sure this is a real match (see below for P- instead)
% axid = 149;
% cellspecs0{axid} = {8};
% matchids0{axid} = 1101;
% neighborids0{axid} = [1118 2420 976 2417];
% neighborids2_0{axid} = [937 992 618 783 844 1326 1381];

% TODO: Looks okay, but P- below seems closer.
% axid = 183;
% cellspecs0{axid} = {8};
% matchids0{axid} = 5781;
% neighborids0{axid} = [4957 6151 5586];
% neighborids2_0{axid} = [5477 5493 6364 6397 6064 3022 4443];

% TODO: Look into; not sure this is a real match (see above for P+ instead)
% axid = 149;
% cellspecs0{axid} = {2};
% matchids0{axid} = 1100;
% neighborids0{axid} = [1039 857 1141 1341];
% neighborids2_0{axid} = [623 620 617 724 4966 4848 4852 1203 1443 1507 993 1113 1291];

% TODO: Looks good, but correlations not nice
axid = 164;
cellspecs0{axid} = {2};
matchids0{axid} = 5176;
neighborids0{axid} = [5056 5317 5316 5362 3951];
neighborids2_0{axid} = [4898 4786 5392 5597 5647 5408];

% TODO: Trace more?  Looks good though.
axid = 174;
cellspecs0{axid} = {2};
matchids0{axid} = 5392;
neighborids0{axid} = [4609 5317 5316];
neighborids2_0{axid} = [5597 5362 5176 5056 4966 486 5206 5628 5507 6486 6366 6121];

% TODO: Looks good, corr needs work (slightly less good P+ match above)
% 4536 Looks contaminated
axid = 183;
cellspecs0{axid} = {2};
matchids0{axid} = 5764;
neighborids0{axid} = [2601];
neighborids2_0{axid} = [5408 5647 5930];


for axid = 1:length(matchids0)
    if exist('eicorrmatch0', 'var') && length(eicorrmatch0) >= axid && ~isempty(eicorrmatch0{axid})
        continue
    end
    
    m = matchids0{axid};
    if isempty(m), continue; end
    axid
    
    n = neighborids0{axid};
    n2 = neighborids2_0{axid};
    cs = cellspecs0{axid};
    eicorrmatch0{axid}   = qd_compare_ei_fit(datarun, axons, axid, matchids0{axid}, 'verbose', true);
    eicorrnear0{axid}    = qd_compare_ei_fit(datarun, axons, axid, [neighborids0{axid} neighborids2_0{axid}], 'verbose', true);
    eicorrall0{axid}     = qd_compare_ei_fit(datarun, axons, axid, cs, 'verbose', true);
    eicorrnotnear0{axid} = setdiff(eicorrall0{axid}, [eicorrmatch0{axid}; eicorrnear0{axid}]);
end


% Add to existing plot from 2007-09-18-4
axonids0 = [174    183  164];
xs0      = [  5    5.5    6];

for i = 1:length(axonids0)
    axid = axonids0(i);
    x = xs0(i);
    plot(randn(size(eicorrnotnear0{axid}))/30 + x, eicorrnotnear0{axid}, '.', 'MarkerSize', 3, 'Color', 0.5*[1 1 1]);
    plot(x, eicorrnear0{axid}, 'ko', 'MarkerSize', 4);
    h = plot(x, eicorrmatch0{axid}, '.', 'MarkerSize', 15, 'Color', [1 0.25 0.25]);
end

textheight = 0.98;
text(xs0(2), textheight, 'OFF Parasol', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
plot(xs0([1 3]), (textheight-0.02)*[1 1], 'k--', 'Clipping', 'off');


%% Find electrodes to cell-finder, hopefully find more cells
cids = [1098 1101 5133 5176 5270];
elecs = arrayfun(...
    @(cid) edu.ucsc.neurobiology.vision.io.NeuronFile.getNeuronIDElectrode(cid), ...
    cids);

radius = 2;
adjelecs = arrayfun(...
    @(e) ai.electrodeMap.getAdjacentsTo(e, radius), ...
    elecs, 'UniformOutput', false);

adjelecs = unique(cell2mat(adjelecs(:)));