


startup_analyse_tenessee
%%


dataset='2014-11-24-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000/data000',dataset);
imov=4;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);

vision_id=1234;

no_cells=5;
noise=0.25;

run_ss

clear all


%%
dataset='2014-11-24-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000/data000',dataset);
imov=4;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);

vision_id=1234;

no_cells=5;
noise=0.25;
 

save_file= sprintf('/Volumes/Lab/Users/bhaishahster/SS_dataset_%s_cell%d_data00%d_a_%dcells.mat',dataset,vision_id,imov,no_cells);

fcn_run_ss(analysis_datafile,bin_datafile,vision_id,no_cells,noise,save_file)

clear all
%%

dataset='2014-11-24-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000/data000',dataset);
imov=4;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);

vision_id=2314;

no_cells=5;
noise=0.25;

run_ss

clear all

%%

dataset='2014-11-24-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000/data000',dataset);
imov=4;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);

vision_id=1039;

no_cells=5;
noise=0.25;

run_ss

clear all

%%

dataset='2014-11-24-0';

%analysis_datafile = sprintf('%s/data001-from-data000/data001-from-data000',dataset);
analysis_datafile = sprintf('%s/data000/data000',dataset);
imov=4;
bin_datafile=sprintf('/Volumes/Data/%s/data00%d',dataset,imov);

vision_id=1621;

no_cells=5;
noise=0.25;

run_ss

clear all