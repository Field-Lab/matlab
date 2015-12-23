%% Not many cones, but not terrible

% Could improve x-scaling
piece = pieces('2012-09-06-0');
conerun = piece.d04s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(13);
rasterrun = piece.d08;
stablerun = piece.d09;
rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 4;
ar = 1.4;
badregions = [];


piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(3);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 3.2;
ar = 1.4;
badregions = [14 16 20];


piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(13);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 4;
ar = 1;
badregions = [15:20];


piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(18);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 5;
ar = 0.9;
badregions = [14 19 20];


piece = pieces('2012-08-21-2');
conerun = piece.d01s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(19)
rasterrun = piece.d03;
stablerun = piece.d04;
rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;
ar = 1;
badregions = [];


piece = pieces('2012-08-21-2');
conerun = piece.d01s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(30)
rasterrun = piece.d03;
stablerun = piece.d04;
rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;
ar = 1;
badregions = [];


piece = pieces('2012-08-21-2');
conerun = piece.d01s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(95)
rasterrun = piece.d03;
stablerun = piece.d04;
rasterrun.rgcs = rasterrun.mapd01s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd01s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 4;
ar = 1;
badregions = [];


%% Pretty bad...
piece = pieces('2012-04-13-1');
conerun = piece.d06s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(20)
rasterrun = piece.d09cm;
stablerun = piece.d10cm;
stablerun.rgcs = stablerun.mapd06s(get_cell_indices(conerun, conerun.rgcs));
rasterrun.rgcs = stablerun.rgcs;
padfactor = 4;
ar = 1.4;
badregions = [];


%% Xscaling didn't work
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(7);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 5;
ar = 1;
badregions = [];


piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(24);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 5;
ar = 1;
badregions = [];


piece = pieces('2012-04-13-1');
conerun = piece.d06s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(13)
rasterrun = piece.d09cm;
stablerun = piece.d10cm;
stablerun.rgcs = stablerun.mapd06s(get_cell_indices(conerun, conerun.rgcs));
rasterrun.rgcs = stablerun.rgcs;
padfactor = 6;
ar = 1;
badregions = [];


piece = pieces('2012-09-06-0');
conerun = piece.d04s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(5);
rasterrun = piece.d08;
stablerun = piece.d09;
rasterrun.rgcs = rasterrun.mapd04s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd04s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 3.8;
ar = 1.3;
badregions = [];


%% Seems like good regions, but terrible RF, xscaling doesn't work
% piece = pieces('2011-12-13-2');
% conerun = piece.d08s;
% conerun.rgcs = conerun.cell_types{2}.cell_ids(4);

%% Seems like good regions, but not good responses.  

% Should cut some spurious cones
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(5);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 5;
ar = 0.75;
badregions = [];


piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(16);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;
ar = .7;
badregions = [];


piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(28);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 2.8;
ar = .7;
badregions = [10];


piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(31);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 3;
ar = 1;
badregions = [];


%% Doubles
piece = pieces('2012-09-21-2');
conerun = piece.d09;
conerun.rgcs = conerun.cell_types{2}.cell_ids(20);
rasterrun = piece.d13;
stablerun = piece.d14;
rasterrun.rgcs = rasterrun.mapd09(get_cell_indices(conerun, conerun.rgcs));
stablerun.rgcs = stablerun.mapd09(get_cell_indices(conerun, conerun.rgcs));
padfactor = 6;
ar = 1;
badregions = [];

%% Not good enough
piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(10);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 3;
ar = .9;
badregions = [];


piece = pieces('2011-12-13-2');
conerun = piece.d08s;
conerun.rgcs = conerun.cell_types{2}.cell_ids(12);
rasterrun = piece.d10;
stablerun = piece.d14s;
rasterrun.rgcs = rasterrun.mapd08s(get_cell_indices(conerun, conerun.rgcs)); rasterrun.rgcs
stablerun.rgcs = stablerun.mapd08s(get_cell_indices(conerun, conerun.rgcs));
padfactor = 4.5;
ar = 1;
badregions = [15];