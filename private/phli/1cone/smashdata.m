%% Load data
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);

d.d01 = load_data('2012-09-24-5/data001', loadopts);
d.d03 = load_data('2012-09-24-5/data003', loadopts);
% d.d05 = load_data('/snle/acquisition/2012-09-24-1/data005/data005', loadopts);
% d.d07 = load_data('/snle/acquisition/2012-09-24-1/data007/data007', loadopts);
% d.d10 = load_data('/snle/acquisition/2012-09-24-1/data010/data010', loadopts);

cellspec = {1,2,3,4,5};

for runname = fieldnames(d)'
    run = d.(runname{1});
    run = load_sta(run, 'load_sta', cellspec, 'guess_stimulus', false);
    run = set_polarities(run);
    d.(runname{1}) = run;
end
clear runname run


%%
smush.names.nickname = 'smushed';
smush.names.short_name = 'smushed';
dnames = fieldnames(d);
firstd = d.(dnames{1});

% Load common fields; have to be the same!
smush.piece = firstd.piece;

% This has to be mostly the same
smush.stimulus = firstd.stimulus;

smush.cell_ids = [];
smush.cell_types = {};
smush.stas.stas = {};
smush.stas.rfs = [];
smush.vision.sta_fits = {};
for i = 1:length(cellspec)
    smush.cell_types{cellspec{i}}.name = firstd.cell_types{cellspec{i}}.name;
    smush.cell_types{cellspec{i}}.cell_ids = [];
end

maxid = 0; % Use to artificially bump up cellids to avoid duplicates
for runname = fieldnames(d)'
    run = d.(runname{1});
    
    smush.cell_ids        = [smush.cell_ids         run.cell_ids+maxid];
    smush.stas.stas       = [smush.stas.stas;       run.stas.stas];
    smush.stas.rfs        = [smush.stas.rfs;        run.stas.rfs];
    smush.vision.sta_fits = [smush.vision.sta_fits; run.vision.sta_fits];
    
    for i = 1:length(cellspec)
        typei = cellspec{i};
        smush.cell_types{typei}.cell_ids = [smush.cell_types{typei}.cell_ids run.cell_types{typei}.cell_ids+maxid];
    end
    clear i typei
    
    maxid = max(smush.cell_ids);
end
clear runname run maxid


%%
datarun = smush;
leave datarun cellspec