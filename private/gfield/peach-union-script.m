clear

data_num = 10;


sub_path{1} = 'done-1/primate-3-5-40-3/data001';
sub_path{2} = 'done-1/primate-3-5-40-4/data001';
sub_path{3} = 'done-1/primate-5-5-40-3/data001';
sub_path{4} = 'done-1/primate-5-5-40-4/data001';
sub_path{5} = 'done-2/primate-3-4-60-3/data001';
sub_path{6} = 'done-2/primate-4-4-60-2-2-15/data001';
sub_path{7} = 'done-2/primate-4-5-40-3-2-15/data001';
sub_path{8} = 'done-2/primate-5-4-60-2-2-15/data001';
sub_path{9} = 'done-2/primate-5-4-60-3/data001';

sub_path{10} = 'data001/data001';
datarun = cell(10,1);

for dset = 1:data_num
    dset
    temp_path = ['/Analysis/gfield/2008-08-27-0/data001/',sub_path{dset}];
    datarun{dset} = load_data(temp_path);
    datarun{dset} = load_sta(datarun{dset}, 'load_sta', []);
    datarun{dset} = load_params(datarun{dset});
    datarun{dset} = get_sta_fits_from_vision(datarun{dset}, 'all');
    
    plot_rf_summaries(datarun{dset}, {1}, 'foa', 1, 'clear', false, 'plot_fits', true')
    plot_rf_summaries(datarun{dset}, {2}, 'foa', 2, 'clear', false, 'plot_fits', true')
    plot_rf_summaries(datarun{dset}, {3}, 'foa', 3, 'clear', false, 'plot_fits', true')
    plot_rf_summaries(datarun{dset}, {4}, 'foa', 4, 'clear', false, 'plot_fits', true')
    plot_rf_summaries(datarun{dset}, {5}, 'foa', 5, 'clear', false, 'plot_fits', true')

end

% compare to original off midget mosaic

orig_datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001');
orig_datarun = load_sta(orig_datarun, 'load_sta', []);
orig_datarun = load_params(orig_datarun);
orig_datarun = get_sta_fits_from_vision(orig_datarun, 'all');

plot_rf_summaries(orig_datarun, {1}, 'foa', 1, 'clear', false, 'plot_fits', true','fit_color', 'r')
title('ON parasol')
plot_rf_summaries(orig_datarun, {2}, 'foa', 2, 'clear', false, 'plot_fits', true','fit_color', 'r')
title('OFF parasol')
plot_rf_summaries(orig_datarun, {3}, 'foa', 3, 'clear', false, 'plot_fits', true','fit_color', 'r')
title('ON midget')
plot_rf_summaries(orig_datarun, {4}, 'foa', 4, 'clear', false, 'plot_fits', true','fit_color', 'r')
title('OFF midget')
plot_rf_summaries(orig_datarun, {5}, 'foa', 5, 'clear', false, 'plot_fits', true','fit_color', 'r')
title('small bistratified')
