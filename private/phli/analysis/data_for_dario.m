%% 120 um
datarun = load_data('2009-12-03-0/data003-nwpca/daat003/daat003');
datarun = load_params(datarun);
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200912030_data003.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200912030_data003.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200912030_data003.ON_Parasol)
plot_ellipses(piece200912030_data003.OFF_Parasol)

clear datarun ans
save p200912030_d003


%% 120 um
datarun = load_data('2012-07-26-0/data000-nwpca');
datarun = load_params(datarun);
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece201207260_data000.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece201207260_data000.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece201207260_data000.ON_Parasol)
plot_ellipses(piece201207260_data000.OFF_Parasol)

clear datarun ans
save p201207260_d000

%%
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2006-06-02-1/data003/data003')
datarun = load_params(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200606021_data003.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200606021_data003.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200606021_data003.ON_Parasol)
plot_ellipses(piece200606021_data003.OFF_Parasol)

clear datarun ans
save p200606021_d003

%% Check same piece as above with bigger stixels
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2006-06-02-1/data000-nwpca/data000/data000')
datarun = load_params(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200606021_data000.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200606021_data000.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200606021_data000.ON_Parasol)
plot_ellipses(piece200606021_data000.OFF_Parasol)

clear datarun ans
save p200606021_d000



%%
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2005-08-03-0/data001/data001')
datarun = load_params(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200508030_data001.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200508030_data001.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200508030_data001.ON_Parasol)
plot_ellipses(piece200508030_data001.OFF_Parasol)

clear datarun ans
save p200508030_d001

%% Check same piece as above with bigger stixels
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2005-08-03-0/data008-nwpca/data008/data008')
datarun = load_params(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200508030_data008.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200508030_data008.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200508030_data008.ON_Parasol)
plot_ellipses(piece200508030_data008.OFF_Parasol)

clear datarun ans
save p200508030_d008

% This run has one OFF parasol that above doesn't have, so maybe worth
% merging that one in.  Cell ID 5178, cell number 455, OFF parasol number
% 49.  In order to combine, the coordinates for this run have to be
% multiplied by 2 as the coordinates are in stixels and the ratio between
% this run and the above is 20:10.



%%
datarun = load_data('2007-02-06-6/data007-nwpca/data007/data007')
datarun = load_params(datarun)
datarun.cell_types = load_txt_cell_types('/snle/lab/Experiments/Array/Analysis/2007-02-06-6/data007-nwpca/data007/data007-parasol-classification.txt', 'order_cell_types', false)
plot_rf_fit(datarun, {2})
plot_rf_fit(datarun, {4})

piece200702066_data007.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
piece200702066_data007.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {4}))
plot_ellipses(piece200702066_data007.ON_Parasol)
plot_ellipses(piece200702066_data007.OFF_Parasol)

clear datarun ans
save p200702066_d007

%% Same piece as above with bigger stixels
datarun = load_data('2007-02-06-6/data005-gdf/data005')
datarun = load_params(datarun)
datarun.cell_types = load_txt_cell_types('/snle/lab/Experiments/Array/Analysis/2007-02-06-6/data005-gdf/data005-parasol-classification.txt', 'order_cell_types', false)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {5})

piece200702066_data005.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200702066_data005.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {5}))
plot_ellipses(piece200702066_data005.ON_Parasol)
plot_ellipses(piece200702066_data005.OFF_Parasol, 'r')

clear datarun ans
save p200702066_d005


%%
datarun = load_data('2008-04-30-1/data005/data005')
datarun = load_params(datarun)
datarun.cell_types = load_txt_cell_types('/snle/lab/Experiments/Array/Analysis/2008-04-30-1/data005/data005-parasol-classification.txt', 'order_cell_types', false)
info(datarun)
plot_rf_fit(datarun, {4})
plot_rf_fit(datarun, {2})

piece200804301_data005.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {4}))
piece200804301_data005.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200804301_data005.ON_Parasol)
plot_ellipses(piece200804301_data005.OFF_Parasol, 'r')

clear datarun ans
save p200804301_d005


%% Parasols
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2006-06-15-0/data000/data000')
datarun = load_params(datarun); info(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece200606150_data000.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece200606150_data000.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece200606150_data000.ON_Parasol)
plot_ellipses(piece200606150_data000.OFF_Parasol, 'r')

clear datarun ans
save p200606150_d000

% ON Midgets
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2006-06-15-0/data002/data002');
datarun = load_params(datarun); info(datarun)
plot_rf_fit(datarun, {3})

piece200606150_data002.ON_Midget  = datarun.vision.sta_fits(get_cell_indices(datarun, {3}))
plot_ellipses(piece200606150_data002.ON_Midget)

clear datarun ans
save p200606150_d002


%%
datarun = load_data('/jacob/snle/lab/Experiments/Array/Analysis/2006-04-24-2/data000/data000')
datarun = load_params(datarun); info(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})
plot_rf_fit(datarun, {5})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
piece.SBC         = datarun.vision.sta_fits(get_cell_indices(datarun, {5}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')
plot_ellipses(piece.SBC, 'g')

piece200604242_data000 = piece;
clear datarun ans piece
save p200604242_d000


%%
datarun = load_data('2005-04-06-2/data007-nwpca/data007')
datarun = load_params(datarun); info(datarun)
plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200504062_data007 = piece;
clear datarun ans piece
save p200504062_d007


%%
datarun = load_data(['/jacob' server_path '2005-05-31-3/data001/data001'])
datarun = load_params(datarun); 
datarun.cell_types = load_txt_cell_types([datarun.names.rrs_prefix '-parasol-classification.txt']);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200505313_data001 = piece;
clear datarun ans piece
save p200505313_d001


%%
datarun = load_data('2008-11-10-0/data013');
datarun = load_params(datarun);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200811100_data013 = piece;
clear datarun ans piece
save p200811100_d013


%%
datarun = load_data('2010-03-05-2/data000');
datarun = load_params(datarun);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece201003052_data000 = piece;
clear datarun ans piece
save p201003052_d000


%%
datarun = load_data('2006-12-04-0/data000');
datarun = load_params(datarun);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200612040_data000 = piece;
clear datarun ans piece
save p200612040_d000


%%
datarun = load_data('2006-12-04-2/data000');
datarun = load_params(datarun);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200612042_data000 = piece;
clear datarun ans piece
save p200612042_d000


%%
datarun = load_data('2005-05-26-6/data000');
datarun = load_params(datarun);
info(datarun)

plot_rf_fit(datarun, {1})
plot_rf_fit(datarun, {2})

piece.ON_Parasol  = datarun.vision.sta_fits(get_cell_indices(datarun, {1}))
piece.OFF_Parasol = datarun.vision.sta_fits(get_cell_indices(datarun, {2}))
plot_ellipses(piece.ON_Parasol)
plot_ellipses(piece.OFF_Parasol, 'r')

piece200505266_data000 = piece;
clear datarun ans piece
save p200505266_d000


%%
for i = 1:length(piece.OFF_Parasol)
    angle.OFF_Parasol(i) = piece.OFF_Parasol{i}.angle;
end
for i = 1:length(piece.ON_Parasol)
    angle.ON_Parasol(i) = piece.ON_Parasol{i}.angle;
end
for i = 1:length(piece.SBC)
    angle.SBC(i) = piece.SBC{i}.angle;
end


%%
for i=1:length(datarun.cell_ids)
    [as, cc{i}] = ei_axon_speed(datarun, datarun.cell_ids(i), 'threshold', 2.5);
    datarun.axon_speed(i) = as;
    disp([num2str(i) ': ' num2str(as)]);
end