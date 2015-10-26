%% Additivity experiment analysis 


nConditions=7;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';


dataRuns_OFF_additivity = [16,17,19,20,22,24,25];
dataRuns_ON_additivity = [16,18,19,21,23,24,26];
movies_OFF_addivitiy =[1,2,5,6,10,13,14];
movies_ON_additivity = [1,4,5,8,12,13,16];
%% Stimulus difference non triviality

% Load stimuli 


% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies_OFF_addivitiy
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-06-0/Visual/null/%d.rawMovie',imov),rawMovFrames,1);
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
%% 
cellID_select = 5453;
sta_depth=30;

%% Load experiment response data


dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1%8,1;
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
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
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


%% fit ASM
WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';

movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= 12*ones(1,length(cellID_list));
destination= 'pc2015_10_06_0_analysis_fits/1SU'
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS

[stas_t,stas_r] = get_gmlm_sta(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth)
save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');

%% Predict response
cellData = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID_select)])

    convolve=150;
    len = 12000;
    binSz=1/1200;
    pred1=cell(2,1);
    icnt=0;
    clear('ss');
for icond=[1,2,3,4,5];
icnt=icnt+1;
    realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1272);
fr= mean(realResp,1)/(1/120);


movd = condMov{icond}-0.5;
 maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 R2_pl = [];

nTrials=30;
 [pred1{icnt},lam]= predictGMLM_bias_lr(cellData.fitGMLM,maskedMov,nTrials,1);
 pred1{icnt} = pred1{icnt}';
 ss(icnt)=spkCondColl{icond};
end
   h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
    % close all
    %figure;
%     ss = spkCondColl{icond};
   hh=plot_record_prediction(ss,pred1);
% 
%     hh=plot_record_pred12(ss,pred4,pred1);
    

%% Compute non-linearity of cell in an LNP model
sta_depth = 30;
WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 

cellIDs = [6741,2176,6961];
data_nls = compute_nl_for_lnp(WN_datafile,cellIDs,sta_depth,movie_xml,stim_length);

figure;
for icell=1:3
subplot(1,3,icell);
plotyy(data_nls(icell).in/data_nls(icell).inp_sd,data_nls(icell).fr,data_nls(icell).XX/data_nls(icell).inp_sd,data_nls(icell).NN/sum(data_nls(icell).NN));
title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

%% Raster difference metric.

dataRuns = dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035/data014-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
imov=14;
WN_datafile_full = '/Volumes/Analysis/2015-10-06-0/streamed/data014/';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId = 2%8,1;
InterestingCell_vis_id= datarun.cell_types{cellTypeId}.cell_ids; %[556,1278,1384,407,1516,2150,2401,3361,4066,4683,5611,6106,6005,7246,7562,3946];
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

spkCondColl=cell(length(dataRuns),1);

condDuration=10;
nConditions=1;
compareConds = [3,4];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-06-0/d14_35-norefit/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%d-from-data014_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
spks = cell(nConditions,1);
for idata = compareConds
spks{idata} = makeSpikeMat(spkCondColl{idata}.spksColl,1/1000,10000);
end

convolve=100;
dist(ref_cell_number) = compare_rasters(spks{compareConds(1)}, spks{compareConds(2)},1/1000,convolve)

end

figure;
hist(dist);

