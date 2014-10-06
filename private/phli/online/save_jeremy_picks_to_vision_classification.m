%% Write Jeremy's picks to vision classification
load_opts = struct('load_params', true);
datarun = load_data('/snle/acquisition/2012-08-21-0/data001/data001', load_opts);

offM = [31
    376
    513
    531
    664
    693
    768
    1009
    1788
    1996
    2973
    3153
    4399
    4831
    5087
    5431
    6953
    6724
    4703];

all = setdiff(datarun.cell_ids, [offM]);

classification = struct('name', 'All', 'cells', all);
classification.subclasses(1) = struct('name', 'Off Midgets', 'cells', offM, 'subclasses', struct([]));
write_txt_classification(classification, [datarun.names.rrs_prefix '_freeman_classification.txt']);