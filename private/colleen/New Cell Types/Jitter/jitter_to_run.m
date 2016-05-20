
clear
close all
%% 2016-04-21-8 data022-cf/edited/data022-cf
dataparam.date='2016-04-21-8/data022-cf/edited';
dataparam.concatname='data022-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-40-1-0.48-11111-20x15-60.35.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[3077 4337 4761 1711 3381 4248 406 3721 407 5912 3140 154 3259 4565 901 3813 7501 1953 4487];
sta = run_jitter(dataparam);

clear 
close all
%% 2016-04-21-8 data022-mVision/data022recluster/ dsata022
dataparam.date='2016-04-21-8/data022-mVision/data022_recluster';
dataparam.concatname='data022';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-40-1-0.48-11111-20x15-60.35.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[4307 4757 5881 1716 422 427 429 4818 4820 4826 4827 152 768 3263 3264];
sta = run_jitter(dataparam);


clear 
close all

%% 2016-04-21-8 data022-mVision/data022_regroup422_4818 data022

dataparam.date='2016-04-21-8/data022-mVision/data022_regroup422_4818';
dataparam.concatname='data022';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-40-1-0.48-11111-20x15-60.35.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[422 770 4827];
sta = run_jitter(dataparam);

clear 
close all
%% 2016-04-21-1 data005-cf/edited/data005-cf

dataparam.date='2016-04-21-1/data005-cf/edited';
dataparam.concatname='data005-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-2-0.48-22222-119.5.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[993 1955 2986 4201 5194 6391 6771 6968 7595 1279 2403 2581 3380 5104 5211 6128 6562 7326 7432 3124 4609 4684 4866 7023 7055 7221];
sta = run_jitter(dataparam);
clear 
close all

%% 2016-04-21-1 data006-cf/edited/data006-cf

dataparam.date='2016-04-21-1/data006-cf/edited';
dataparam.concatname='data006-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-22222-119.5.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[527 902 1622 4595 5193 5493];
sta = run_jitter(dataparam);
clear 
close all

%% 2015-09-23-7 data031-cf/edited/data031-cf

dataparam.date='2015-09-23-7/data031-cf/edited';
dataparam.concatname='data031-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-8-0.48-11111.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[526 1022 1568 2371 3186 4041 4703 5239 6781 7101 917 3931 5386 6827 7216 306 413 996 1051 1759 1816 3396 3586 3830 4656 5733 6061 7148 81 800 966 1174 2912 3063 3187 3376 3683 3768 4088 4861 5031 5272 5418 6199 6323 7127 7563];
sta = run_jitter(dataparam);


clear 
close all

%% 2015-09-23-7 data028-cf/edited/data028-cf

dataparam.date='2015-09-23-7/data028-cf/edited';
dataparam.concatname='data028-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-20-12-0.48-11111.xml';
dataparam.select_cells = 1;
dataparam.cell_specification =[526 527 528 1027 1430 1568 2377 3171 3215 4201 4202 4203 4714 5238 6786 6967 875 966 3931 4342 4703 6368 7328 305 412 1876 1006 1007 1008 3498 4651 5341 5660 5941 7043 151 968 2915 3576 4792 5227 5418];
sta = run_jitter(dataparam);

