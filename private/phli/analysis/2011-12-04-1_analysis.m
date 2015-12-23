%% Basics
piece = '2011-12-04-1';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
udprintpath = [printpath() '/ud'];
allconesprintpath = [printpath() '/allcones'];
keep_vars = {'piece'; 'loadopts'; 'udprintpath'; 'allconesprintpath'; 'keep_vars'; 'd'};


%% Streamed runs and cell picks
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};

d.d00s = load_data([piece '/streamed/data000-0/data000-0'], streamedopts);
d.d00s.offP = [61 293 1237 1561 7007];

d.d02s = load_data([piece '/streamed/data002-0/data002-0'], streamedopts);
d.d02s.offP = [34 289 1636];%[34 289 1158 1636 7157];


%% data004 C/R 
if ~isfield(d, 'd04'), d.d04 = load_data([piece '/data004'], loadopts); end
d.d04 = read_stim_lisp_output(d.d04, '2011-12-04-1_f00_1234serial');
d.d04.stimulus = convert_stimulus_to_combined_maps(d.d04.stimulus);
d.d04.stimulus = parse_stim_rgbs(d.d04.stimulus);
d.d04.mapd02s = map_ei(d.d02s, d.d04);

conerun = d.d00s;
conerun.rgcs = conerun.offP;

% Since d00s didn't have EIs and I can't pull data from rush today, hack
% around using matching cells from d02s.
rasterrun = d.d04;
rasterrun.rgcs = rasterrun.mapd02s(get_cell_indices(d.d02s, [34 289 1158 1636 7157]));

[rasterrun, conerun, ~, crs, crsx, ps, resnorms, residuals] = cr_compound_plot(rasterrun, conerun, rasterrun.triggers(1:2:end), 'rgcindices', 'all');


%% Stability summary
piece = '2011-12-04-1';

d0 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, 0, 0));
d0 = load_globals(d0);
d0time = d0.globals.readRDH512().getTime();

coneruns = [0 2 6 8];
numplots = length(coneruns) - 1;

fig = figure();
for i = 1:numplots
    subplot(1, 3, i);
    
    d1run = coneruns(i);
    d2run = coneruns(i+1);

    d1 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, d1run, d1run));
    d1 = load_params(d1);
    d1 = load_cones(d1);
    d1 = load_globals(d1);
    d1 = stimulus_from_globals(d1);
    
    d2 = load_data(sprintf('%s/streamed/data%.3d-0/data%.3d-0', piece, d2run, d2run));
    d2 = load_params(d2);
    d2 = load_cones(d2);
    d2 = load_globals(d2);
    d2 = stimulus_from_globals(d2);
    
    overlay_cone_mosaics(d1, d2);
    axis equal tight
    
    secdiff = d2.globals.readRDH512.getTime - d1.globals.readRDH512.getTime;
    startdiff = d2.globals.readRDH512.getTime - d0time;
    title(sprintf('data%.3d vs data%.3d\n%s apart\n%s since data000', d1run, d2run, printtime(secdiff), printtime(startdiff)));
end


% Print
set(fig, 'PaperType', 'B');
orient(fig, 'landscape');
set(fig, 'Color', 'w');
set(fig,'InvertHardcopy','off')
print('-dpdf', '-painters', '-noui', printpath('stability', sprintf('stability-%s', piece)));