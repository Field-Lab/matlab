


% Nishal P. Shah , September 2014

% Have a dataset ready with cells you want to nullify grouped under a
% common name(use Vision software).

startup_bertha

startup_rooster

datafile='2014-11-05-9/data001/data001';%'nishal/2014-08-20-2/data001/data001';

% type_name_inp = 'userCellList' for a list of cells.

no_images_per_movie=10;
start_image=10;
tag='pc_contrast_analysis';
destination_mat=['/Volumes/Analysis/nishal/',tag];
if(~exist(destination_mat,'dir'))
mkdir(destination_mat);
end

dest_raw='/Volumes/Data/stimuli/movies/null_space/';

movies=cell(20,1);
% 'solver' .. 1 for LSQR, 2 for CRAIG, 3 for Fourier 


%% BW - 10 sec user List

cell_params=struct();
cell_params.type_name_inp='userCellList';%'nc2';%'userCellList';
cell_params.cell_list=[272,5824,5825,6378]%[3888,2825,1820,4129, 5346,5671,5161,1278, 3828,3574,4036,3572, 503,560,797,1009,487,181,901]; % if type_name_inp = 'userCellList' 
cell_params.use_fits=2;
cell_params.STAlen=14;

mov_params=struct();
mov_params.mov_type='bw-precomputed';
mov_params.movie_time=120*10;
mov_params.mean=0.5*255;
mov_params.deviation=0.48*255;
mov_params.scaling_loss=0.05; % a number in [0,1], fraction of values that is changed by scaling.


solver=3;
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params,mov_params,solver);

%%

load('/Volumes/Analysis/nishal/test-null.mat');

mov_len=size(mov_orignial,3);
% Variance for whole movie
sd_orig=[];
sd_null=[];
for icell=1:length(CellMasks)
    mask=logical(repmat(CellMasks{icell},[1,1,mov_len]));
    mov_orig_mask=mov_orignial(mask);
    mov_null_mask=mov_modify_new(mask);
    
    sd_orig=[sd_orig;sqrt(var(mov_orig_mask(:)))];
    sd_null=[sd_null;sqrt(var(mov_null_mask(:)))];
end
% 
% [N1,X1]=hist(sd_orig,20);
% [N2,X2]=hist(sd_null,20);
% figure;
% bar(X1,N1,'r');
% hold on;
% bar(X2,N2,'b');
% legend('Original','Null');
figure('Color','w');
plot(sd_orig,sd_null,'*');
xlabel('Original');
ylabel('Null');
hold on;
plot([0,min(sd_orig),max(sd_orig)],[0,min(sd_orig),max(sd_orig)],'r');
legend('Orig v/s Null','45d line');
title('Overall variance');

% Variance frame wise

% Variance STA- wise