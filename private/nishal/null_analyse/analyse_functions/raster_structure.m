%% ON - OFF assymetry in low contrast nulling??

%% pc 2015-11-09-1

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
a
WN_datafile = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009';
WN_datafile_streamed = '/Volumes/Analysis/2015-11-09-1/streamed/data009/data009';
location = '/Volumes/Analysis/2015-11-09-1/d14_36-norefit';
%% pc2015-10-29-7
% Condition strings
nConditions=6;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [14,15,17,18,20,22,24];
dataRuns_ON_additivity = [14,16,17,19,21,23,24];
movies_OFF_addivitiy =[1,2,5,6,10,14,11];
movies_ON_additivity = [1,4,5,8,12,16,11];


WN_datafile ='/Volumes/Analysis/2015-10-29-7/d12_41-norefit/data012/data012';
WN_datafile_streamed ='/Volumes/Analysis/2015-10-29-7/streamed/data012/data012';
location = '/Volumes/Analysis/2015-10-29-7/d12_41-norefit';

%% pc 2015-10-29-2

% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [3,4,6,7,9,11,13];
dataRuns_ON_additivity = [3,5,6,8,10,12,13];
movies_OFF_addivitiy =[1,2,5,6,10,14,13];
movies_ON_additivity = [1,4,5,8,12,16,13];


WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001';
WN_datafile_streamed = '/Volumes/Analysis/2015-10-29-2/streamed/data001/data001';
location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

%% pc 2015-11-09-8

% Condition strings
nConditions=6;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [15,16,18,19,21,23];
dataRuns_ON_additivity = [15,17,18,20,22,25];
movies_OFF_addivitiy =[1,2,9,6,10,14];
movies_ON_additivity = [1,4,9,8,12,16];


WN_datafile = '/Volumes/Analysis/2015-11-09-8/d15_61-norefit/data007/data007';
WN_datafile_streamed = '/Volumes/Analysis/2015-11-09-8/streamed/data007/data007';
location ='/Volumes/Analysis/2015-11-09-8/d15_61-norefit';

%%


dataRuns =dataRuns_OFF_additivity;
cellTypeId = 2;

datarun_s=load_data(WN_datafile_streamed);
datarun_s=load_params(datarun_s);
datarun_s = load_sta(datarun_s);
InterestingCell_vis_id_streamed=datarun_s.cell_types{cellTypeId}.cell_ids; 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);

InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;
conda=3;condb=4;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
conda_var=[];condb_var=[];conda_mean=[];condb_mean=[];cells_select=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID= InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
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
    
    % compare structure in rasters in conditions a and b
      convolve=150;
    binSz = 1/1200;len=12000;
    realResp = makeSpikeMat(spkCondColl{conda}.spksColl, binSz,len);
    [PSTH_reca,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        
    realResp = makeSpikeMat(spkCondColl{condb}.spksColl,binSz,len);
    [PSTH_recb,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
    
    
    conda_var = [conda_var,sqrt(var(PSTH_reca(0.25*end:end)))];
    condb_var = [condb_var,sqrt(var(PSTH_recb(0.25*end:end)))];
    
    conda_mean = [conda_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
    condb_mean = [condb_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
     
    close all
    
    
   
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(6,3,[1:9]);
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    
    subplot(6,3,[10,11,12]);
   plot(time,PSTH_reca);hold on;plot(time,PSTH_recb);
    title(sprintf('%f,%f',conda_var(end),condb_var(end)));

    subplot(6,3,[13,16]);
    imagesc(mean(datarun.stas.stas{datarun.cell_ids==cellID}(:,:,:,27),3));
    colormap gray
    axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
   
    subplot(6,3,[14,17]);
   plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
   hold on;
    plot_rf_fit_nishal(datarun, cellID,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',true,'edge',true,'labels',false,'edge_color',[1,0,0]);
   axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
  title('Concat class.');
      
    subplot(6,3,[15,18]);
   plot_rf_fit_nishal(datarun_s, InterestingCell_vis_id_streamed,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
   axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
 title('Streamed class.');
 
   suptitle(sprintf('CellID: %d',cellID));
% dialogue if a good cell or a bad cell
% Construct a questdlg with two options
cell_ch=1;
choice = questdlg('Is it a good cell?', ...
	'Cell choice', ...
	'good cell','bad cell','good cell');
% Handle response
switch choice
    case 'good cell'
        disp([choice ' selected'])
        cell_ch = 1;
    case 'bad cell'
        disp([choice ' selected'])
        cell_ch = 0;
end


cells_select = [cells_select;cell_ch];
end

condba_rat = condb_var./conda_var;

data(cellTypeId).conda_var=conda_var;
data(cellTypeId).condb_var=condb_var;
data(cellTypeId).condba_rat=condba_rat;
data(cellTypeId).cells_select = cells_select;
data(cellTypeId).conda_mean =conda_mean;
data(cellTypeId).condb_mean =condb_mean;
data(cellTypeId).cellIDs = InterestingCell_vis_id;

figure;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 
 %% plot histogram
 load('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_8/d_add/PSTH_cond3_4.mat');
 figure;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 title('pc 2015 11 09 8')
 
  load('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_add/PSTH_cond3_4.mat');
 figure;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 title('pc 2015 11 09 1')
 
load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
 figure;
 histogram(data(1).condba_rat(data(1).cells_select==    1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 title('pc 2015 10 29 2')
 
  load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_7/d_add/PSTH_cond3_4.mat');
 figure;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 title('pc 2015 10 29 7')
 
 %%
  load('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_8/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.45; lb= 0.4; 
  
cls = data(cellType).cellIDs(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 
cols='rkrkrkrkrkrkkrkrkrkr';
cellID=cls(6)
       
condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
    h=figure('units','normalized','Position',[0 0 0.2 0.125],'Color','w');
    for irun = 3:4%1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    xlim([2.5,10]);
    set(gca,'yTick',[]);axis off
    %xlabel('Time (s)->')
    
    
    
    %% plot non-linearities for pc 2015-10-29-2
      load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.25; lb= 0.2; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data014/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data014/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
    figure;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'b');
h1.Color(4)=0.5;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.5; lb= 0.4; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data014/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data014/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
hold on;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'r');
h1.Color(4)=0.3;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

%% conditions WN low contrast and WN+Null


dataRuns = dataRuns_OFF_additivity;
cellTypeId=1;
 
  data34 = load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
  

datarun_s=load_data(WN_datafile_streamed);
datarun_s=load_params(datarun_s);
datarun_s = load_sta(datarun_s);
InterestingCell_vis_id_streamed=datarun_s.cell_types{cellTypeId}.cell_ids; 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);

InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;
conda=3;condb=5;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
conda_var=[];condb_var=[];conda_mean=[];condb_mean=[];cells_select=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID= InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    
    % compare structure in rasters in conditions a and b
      convolve=150;
    binSz = 1/1200;len=12000;
    realResp = makeSpikeMat(spkCondColl{conda}.spksColl, binSz,len);
    [PSTH_reca,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        
    realResp = makeSpikeMat(spkCondColl{condb}.spksColl,binSz,len);
    [PSTH_recb,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
    
    
    conda_var = [conda_var,sqrt(var(PSTH_reca(0.25*end:end)))];
    condb_var = [condb_var,sqrt(var(PSTH_recb(0.25*end:end)))];
    
    conda_mean = [conda_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
    condb_mean = [condb_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
     
    

end

condba_rat = condb_var./conda_var;

data(cellTypeId).conda_var=conda_var;
data(cellTypeId).condb_var=condb_var;
data(cellTypeId).condba_rat=condba_rat;
data(cellTypeId).cells_select = data34.data(cellTypeId).cells_select;
data(cellTypeId).conda_mean =conda_mean;
data(cellTypeId).condb_mean =condb_mean;
data(cellTypeId).cellIDs = InterestingCell_vis_id;

figure;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
 	save('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_5.mat','data','conda','condb');
