addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha


%% Additivity experiment analysis 


% Condition strings
nConditions=6;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [14,15,17,18,20,22];
dataRuns_ON_additivity = [14,16,17,19,21,23];
movies_OFF_addivitiy =[1,2,9,6,10,14];
movies_ON_additivity = [1,4,9,8,12,16];

%% Stimulus difference non triviality

% Load stimuli 


% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies_OFF_addivitiy
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
end

% Display stimulus
figure;
for itime=1:100
subplot(3,1,1);
    imagesc(condMov{3}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot(3,1,2);
    imagesc(condMov{4}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot (3,1,3);
    imagesc(abs(condMov{3}(:,:,itime)' - condMov{4}(:,:,itime)'));axis image;colormap gray;colorbar;caxis([0,1]);
pause(1/120);
end

% Display stimuli with cell locations

% Euclidean difference of movies
%% Predict response  and compare with real response! 

% select cell
cellID = 5341%cellID_select(1);

% set up stuff for reading real responses
dataRuns = dataRuns_ON_additivity;
WN_datafile ='/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009'

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end

% predict response
%%

% destination = 'pc2015_11_09_1_analysis_fits/SUs_data025';
% cellData001 = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)]);

nSU=5;        
destination = 'pc2015_11_09_1_analysis_fits/SUs_data025_2';
cellData = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)])
fitGMLM = cellData.fitGMLM_log{nSU};
    
scale_tf=1;

pred1=cell(2,1);
clear('ss');
icnt=0;
for icond=[1,2]%[1,2,3,4,5];
icnt=icnt+1;
    realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
fr= mean(realResp,1)/(1/120);


movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
movd_masked = movd(repmat(cellData.mask>0,[1,1,size(movd,3)]));
movd=movd/3;%*0.4804/sqrt(var(movd_masked(:)));
 maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 R2_pl = [];

nTrials=100;
 [pred1{icnt},lam]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,1);
%  [pred1{icnt},lam]= predictGMLM_gamma2(fitGMLM,maskedMov2,nTrials,2,1);
 pred1{icnt} = pred1{icnt}';
 ss(icnt)=spkCondColl{icond};
end
%close all;

h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
hh=plot_record_prediction(ss,pred1);
xlim([8,10]);

%% Load experiment response data


dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 8%8,1;
InterestingCell_vis_id=cellID_select% datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_10_06_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
%     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s',datarun.cell_types{cellTypeId}.name)))
%         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s',datarun.cell_types{cellTypeId}.name));
%     end
%    s=hgexport('readstyle','raster');
%    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_06_0/dsu_null_16,17,19,20,22,24,25/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

%% Difference of stimulus in additivity experiment..

dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun)

cellTypeId = 2%8,1;
InterestingCell_vis_id= datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cond1 =(3);cond2=(4);

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
normmov1=[];normmov2=[];normmov_diff=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
  close all
  cellID=InterestingCell_vis_id(ref_cell_number)
  stas=datarun.stas.stas(datarun.cell_ids==cellID);
  
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:size(stas{1},4)
         st_temp(:,:,:,itime)=mean(stas{1}(:,:,:,end-itime+1),3)'; 
    end
    %sprintf('Only blue gun selected!!!!!')
    stas{1}=st_temp;
 
cell_params.STAlen=14;
  [new_stas,totalMaskAccept,CellMasks]=clipSTAs(stas,cell_params)

  mov1 = (condMov{cond1}-0.5).*repmat(totalMaskAccept,[1,1,size(condMov{cond1},3)]);
   mov2 = (condMov{cond2}-0.5).*repmat(totalMaskAccept,[1,1,size(condMov{cond2},3)]);
    mov1 = mov1-mov2;
    
    normmov1(ref_cell_number) = norm(mov1(:));
    normmov2(ref_cell_number) = norm(mov2(:));
    normmov_diff(ref_cell_number) = norm(mov1(:)-mov2(:));
  
end

figure;
plot(normmov1,normmov_diff,'.');
hold on;
plot([0,max(normmov_diff)],[0,max(normmov_diff)],'g');
hold on;
plot([0,max(normmov_diff)],[0,2*max(normmov_diff)],'g');
xlabel('low contrast WN');
ylabel('norm of null movie added');

figure;
plot(normmov2,normmov_diff,'.');
hold on;
plot([0,max(normmov_diff)],[0,max(normmov_diff)],'g');
hold on;
plot([0,max(normmov_diff)],[0,2*max(normmov_diff)],'g');
xlabel( 'WN + null');
ylabel('norm of null movie added');


figure;
plot(normmov1,normmov2,'.');
hold on;
plot([0,max(normmov2)],[0,max(normmov2)],'g');
hold on;
plot([0,max(normmov2)],[0,2*max(normmov2)],'g');
xlabel( 'low contrast WN');
ylabel('WN + null');

%% compare WN v/s WN+null stimuli for a particular cell


dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun)

cellTypeId = 2%8,1;
InterestingCell_vis_id= 6741%cellID_select; %datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cond1 =(3);cond2=(4);

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
normmov1=[];normmov2=[];normmov_diff=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
  close all
  cellID=InterestingCell_vis_id(ref_cell_number)
  stas=datarun.stas.stas(datarun.cell_ids==cellID);
  
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:size(stas{1},4)
         st_temp(:,:,:,itime)=mean(stas{1}(:,:,:,end-itime+1),3)'; 
    end
    %sprintf('Only blue gun selected!!!!!')
    stas{1}=st_temp;
 
cell_params.STAlen=14;
  [new_stas,totalMaskAccept,CellMasks]=clipSTAs(stas,cell_params)

  mov1 = (condMov{cond1}-0.5).*repmat(totalMaskAccept,[1,1,size(condMov{cond1},3)]);
   mov2 = (condMov{cond2}-0.5).*repmat(totalMaskAccept,[1,1,size(condMov{cond2},3)]);
    mov1 = mov1-mov2;
    
    normmov1(ref_cell_number) = norm(mov1(:));
    normmov2(ref_cell_number) = norm(mov2(:));
    normmov_diff(ref_cell_number) = norm(mov1(:)-mov2(:));
  
end
[r,c] = find(totalMaskAccept==1)

figure;
for itime=100:200%1:size(mov1,3)
subplot(1,3,1);
imagesc(mov1(min(r):max(r),min(c):max(c),itime)');
colormap gray;
axis image
caxis([-0.5,0.5]);
set(gca,'xTick',[]); set(gca,'yTick',[]);
title('low contrast WN');

subplot(1,3,2);
imagesc(mov2(min(r):max(r),min(c):max(c),itime)');
colormap gray;
axis image
caxis([-0.5,0.5]);
set(gca,'xTick',[]); set(gca,'yTick',[]);
title('WN+ NUll movie')

subplot(1,3,3);
imagesc(mov2(min(r):max(r),min(c):max(c),itime)' - mov1(min(r):max(r),min(c):max(c),itime)');
colormap gray;
axis image
caxis([-0.5,0.5]);
set(gca,'xTick',[]); set(gca,'yTick',[]);
title('Null movie added');
pause(1/50)
end



%% Analyse a cell - fit ASM, difference of rasters, non-linearitiy ..

%% 
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 1%OFF:2 .. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]


 %% fit ASM
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data026/data026'

movie_xml = 'BW-8-4-0.48-11111';
stim_length=900;% 
cellID_list = [2477,511,5867]%[cellID_select];
nSU_list= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
destination= 'pc2015_11_09_1_analysis_fits/SUs_data026_exp'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)
save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');



 %% fit ASM
WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001'

movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
cellID_list = 198%[cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22];
destination= 'pc2015_10_29_2_analysis_fits/SUs_data001_quad2'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[~] = get_gmlm_sta3(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth)
save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');




%% Explore time course
icond=2;
    realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
fr= mean(realResp,1)/(1/120);

pred1=cell(2,1);
clear('ss');

icnt=0;

for nSU=15%1:10
    fitGMLM = fitGMLM_log{nSU};%cellData.fitGMLM_log{nSU};
for scale_tf = 0%[1:10:60];%[0.5:0.1:1]
    icnt=icnt+1;
movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
movd_masked = movd(repmat(cellData.mask>0,[1,1,size(movd,3)]));
movd=movd*0.4804/sqrt(var(movd_masked(:)));
 maskedMov= filterMov(movd,cellData.mask,transform_tf(squeeze(cellData.ttf),scale_tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 
nTrials=30;
 %[pred1{icnt},lam]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,1);
 %[pred1{icnt},lam]= predictGMLM_gamma2(fitGMLM,maskedMov2,nTrials,2,1);
  [pred1{icnt},lam]= predictGMLM_gamma_opnl(fitGMLM,maskedMov2,nTrials,1);
  
 pred1{icnt} = pred1{icnt}';
 ss(icnt)=spkCondColl{icond};
end
end

h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
hh=plot_record_prediction(ss,pred1);
%xlim([3,6]);
%%
% have meanR, lam
meanRe= mean(realResp,1);
%% Compute non-linearity of cell in an LNP model
sta_depth = 30;
WN_datafile_cellID ='/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009'
WN_datafile ='/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data026/data026'
movie_xml = 'BW-8-4-0.48-11111';contrast_factor=0.5;
stim_length=900;% 

datarun=load_data(WN_datafile_cellID)
datarun=load_params(datarun)

tic;
cellTypeId = 1%8,1;
destination= 'pc2015_11_09_1_analysis_fits/lnp_data026'
cellIDs = datarun.cell_types{cellTypeId}.cell_ids; % [6741,2176,6961]
data_nls = compute_nl_for_lnp(WN_datafile,cellIDs,sta_depth,movie_xml,stim_length,contrast_factor);
%save(['/Volumes/Lab/Users/bhaishahster/',destination,'/data_nls2.mat'],'data_nls');
toc;

iidx = 1:length(datarun.cell_types{1}.cell_ids);
mcellid=[];
for icell=[2477,511,5867,4726,6437,7636,6978,5566]
mcellid = [mcellid;iidx(datarun.cell_types{1}.cell_ids == icell)];
end

figure;
for icell=mcellid'
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
hold on;
plot(data_nls(icell).in,data_nls(icell).fr/100);
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end


figure;icnt=0;
for icell=mcellid'
    icnt=icnt+1;
subplot(1,3,icnt);
xx=data_nls(icell).NN;
plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,xx/(sum(xx)));
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

% plotNL for different movies 
% have condMov ready
data = data_nls(54);
figure;
for icond=1:6
icell=1;
mov = condMov{icond}-0.5; % make sure it's between -0.5 to 0.5
maskedMovdd= filterMov(mov,data.mask,squeeze(data.ttf));
input = (data.usta'*maskedMovdd - data.ainput) / data.binput;
subplot(3,2,icond);
hist(input,20);
hold on;
plot(data.in,data.fr,'r');
xlim([-10,10]);
end

%% Raster difference metric.

dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 2%8,1;
InterestingCell_vis_id= datarun.cell_types{cellTypeId}.cell_ids; %[6741,2176,6961]
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

spkCondColl=cell(length(dataRuns),1);

condDuration=10;
nConditions=1;
compareConds = [3,6];
dist4_from3=[];dist_sd4_from3=[];
dist3_from4=[];dist_sd3_from4=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end

    
spks = cell(nConditions,1);
for idata = compareConds
spks{idata} = makeSpikeMat(spkCondColl{idata}.spksColl,1/1000,10000);
end 

convolve=100;
[dist4_from3(ref_cell_number),dist_sd4_from3(ref_cell_number)] = compare_rasters(spks{compareConds(2)}, spks{compareConds(1)},1/1000,convolve)
[dist3_from4(ref_cell_number),dist_sd3_from4(ref_cell_number)] = compare_rasters(spks{compareConds(1)}, spks{compareConds(2)},1/1000,convolve)

end

figure;
hist(dist);

%% plot NL indix v/s difference of raster

% plots nl idx for different cells
figure;
for icell=1:length(data_nls)
subplot(6,5,icell);
plot(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).in,data_nls(icell).NL_fit.g(data_nls(icell).in));
title(sprintf('%d:',icell));
end

% plots mask of pixels for different cells
figure;
for icell=1:length(data_nls)
subplot(6,5,icell);
imagesc(data_nls(icell).mask');
colormap gray;
axis image
set(gca,'xTick',[]);set(gca,'yTick',[]);
title(sprintf('%d:',icell));
end

goodCells = 1:length(InterestingCell_vis_id);%[1,2,3,4,6,7,8,9,10,12,13,14,17,20,21,28]
nlidx=[];
dist =[];
for icell=goodCells
nlidx = [nlidx;data_nls(icell).NL_fit.gmed];
dist = [dist;abs(dist3_from4(icell))];
end

figure;
plot(dist,nlidx,'.')
xlabel('Distance- prediction of resp to wn from wn+null');
ylabel('Non-linearity index');

nlidx=[];
dist =[];
for icell=goodCells
nlidx = [nlidx;data_nls(icell).NL_fit.gmed];
dist = [dist;abs(dist4_from3(icell))];
end

figure;
plot(dist,nlidx,'.')
xlabel('Distance- prediction of resp to wn+null from wn');
ylabel('Non-linearity index');

%% make conditional interaction plots
WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002';
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;
cellID=49;
userSTA_depth=30;
destination = 'pc2015_10_29_2_analysis_fits/SUs_data002';
ASM_link =['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)];
nSU=4;

sus=[1,2];
suP=[1,2,3,4];
u1=[1,2];
u2=[3,4];
upossible=[3,4];
[sd,dist,su_ass,sd_self]  = probe_cond_interaction_fcn_pairs(WN_datafile,movie_xml,stim_length,cellID,userSTA_depth,ASM_link,nSU)


%% compute output NL

WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001';
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;
cellID=198;
userSTA_depth=30;
destination = 'pc2015_10_29_2_analysis_fits/SUs_data001_quad2';
ASM_link =['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)];

fitGMLM_log = compute_op_NL_WN(WN_datafile,movie_xml,stim_length,cellID,userSTA_depth,ASM_link)