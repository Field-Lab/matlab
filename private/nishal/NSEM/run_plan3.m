
startup_analyse_tenessee
%%


dataset='2006-07-07-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000-mg/data000/data000',dataset);
imov=0;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);
imov=1;
vision_id=5403;

no_cells=5;
noise=0.25;

run_ss

clear all
