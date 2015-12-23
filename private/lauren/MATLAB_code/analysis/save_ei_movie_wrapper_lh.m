datarun = load_data('2008-08-27-2/data001-lh/data001-lh');

datarun = load_neurons(datarun);
datarun = load_ei(datarun,'all');
datarun.piece.array_id = 60;

save_ei_movie(datarun, 811, '/snle/home/lhruby/Desktop/n811_ei_movie.mov')