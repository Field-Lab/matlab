%% ON - OFF assymetry in low contrast nulling??

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

%% compute structure in raster - code in section 2 of raster_structure..
    
    %% plot non-linearities for pc 2015-10-29-2
    
    %% on v/s off non-linearity
      load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.25; lb= 0.2; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
    figure;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0,0.4470,0.7410],'LineWidth',2);
h1.Color(4)=0.5;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.5; lb= 0.4; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
hold on;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.8500,0.3250,0.0980],'LineWidth',2);
h1.Color(4)=0.2;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

set(gca,'Visible','off')

    %% on v/s off non-linearity - in same range
      load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.4; lb= 0.3; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
    figure;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0,0.4470,0.7410],'LineWidth',2);
h1.Color(4)=0.5;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.4; lb= 0.3; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
hold on;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.8500,0.3250,0.0980],'LineWidth',2);
h1.Color(4)=0.2;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

set(gca,'Visible','off')


%% ON cell type NL analysis
      load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.1; lb= 0.05; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
    figure;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h1.Color(4)=0.5;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.35; lb= 0.3; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
hold on;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.4660,0.6740,0.1880],'LineWidth',2);
h1.Color(4)=0.2;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

set(gca,'Visible','off')
%% Plot OFF NL

      load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.7; lb= 0.6; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
    figure;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h1.Color(4)=0.5;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.4; lb= 0.3; 
 
  if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end
  
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
  
hold on;
for icell=mcellid
%plotyy(data_nls(icell).in,data_nls(icell).fr,data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)));
h1=plot(data_nls(icell).XX,data_nls(icell).NN/(sum(data_nls(icell).NN)),'Color',[0.25,0.25,0.25]);
h1.Color(4)=0.1;
hold on;
h1=plot(data_nls(icell).in,data_nls(icell).fr/100,'Color',[0.4660,0.6740,0.1880],'LineWidth',2);
h1.Color(4)=0.2;
hold on;
%title(sprintf('Cell: %d, nl: %0.02f',cellIDs(icell), data_nls(icell).NL_fit.nl_idx));
%title(sprintf('Cell: %d',cellIDs(icell)));
%xlim([-0.2,0.2]);
end

set(gca,'Visible','off')
%%  plot nlidx  and null response structure 
% select how to parametrize non-linearity
 
load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');


dataRuns = dataRuns_OFF_additivity;
cellType=2;
load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
 
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).cells_select'==1);
  
figure;
for icell=mcellid
    if(data_nls(icell).cellID== data(cellType).cellIDs(icell))
        display('cell IDs match');
    end
    plot(data(cellType).condba_rat(icell),data_nls(icell).NL_fit.g0,'b.');hold on;
end

hold on;

dataRuns = dataRuns_ON_additivity;
cellType=1;
load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
 
  iidx = 1:length(data_nls);
  mcellid = iidx(data(cellType).cells_select'==1);
  

for icell=mcellid
    if(data_nls(icell).cellID==  data(cellType).cellIDs(icell))
        display('cell IDs match');
    end
    plot(data(cellType).condba_rat(icell),data_nls(icell).NL_fit.g0,'r.');hold on;
end

%% prediction using non-linearities
    load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
  cellType=1;ub = 0.5; lb= 0.4; 
 movies = movies_ON_additivity;


% load nls
 if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end


  iidx = 1:length(data_nls);
  mcellid = iidx%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 cids = data(cellType).cellIDs;
  
  
% load movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-2/Visual/%d.rawMovie',imov),rawMovFrames,1);
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

% iterate over different cells
gain_list = [1,2,3,5,7,10,20,30];
pred_log = cell(length(mcellid),length(gain_list),2);
R2_log = zeros(length(mcellid),length(gain_list),2);
cellID_log=[];
 for imcell=1:length(mcellid)
     display(sprintf('Finished: %d out of %d cells',imcell,length(mcellid)));
  current_cell=mcellid(imcell);
  cellID = cids(current_cell);
% load true responses .
cellID_log=[cellID_log;cellID];
   if(data_nls(current_cell).cellID== cellID)
  display('cell IDs match');
   end
cols='rkrkrkrkrkrkkrkrkrkr';

       
condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
    igain_cnt=0;

for gain=gain_list
    igain_cnt=igain_cnt+1;
        R2_pl=[];
    icnt=1;pred=cell(2,1);clear ss;nTrials=30;
    for icond=3:4
        pred{icnt} = response_lnp(data_nls(current_cell),gain*(condMov{icond}-0.5),nTrials);    
        ss(icnt)=spkCondColl{icond};
        
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        convolve=15;
        binSz = 1/120;len=1200;
        [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        pred{icnt}(pred{icnt}>100000)=100000;
        [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});

        R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
        PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
        
        R2_log(imcell,igain_cnt,icnt)=R2_pl(end);
        p{1} =pred{icnt};
        pred_log(imcell,igain_cnt,icnt) =p;
        icnt=icnt+1;
        
    end
    
% h=figure;%('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
% PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
% hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
% title(sprintf('gain: %d, %0.04f,%0.04f',gain,R2_pl(1),R2_pl(2)));
% xlim([8,10]);
end  

end
save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/predictions_cond3_4_ON.mat','pred_log','R2_log','gain_list','cellID_log');
 

%%

load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');


cellType=2;
ub = 0.9; 
lb= 0.4;



iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
mcellid = iidx(cell_bin);
cids = data(cellType).cellIDs(cell_bin);

if(cellType==2)
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/predictions_cond3_4_OFF.mat');
else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/predictions_cond3_4_ON.mat');
end
gain_list = [1,2,3,5,7,10,20,30];

figure;
for iplot=1:8
subplot(2,4,(iplot));
plot(R2_log(cell_bin,iplot,1),R2_log(cell_bin,iplot,2),'.');
xlim([0,1]);ylim([0,1]);
axis square
title(sprintf('gain: %d',gain_list(iplot)));
end


%% plot recorded and predicted responses

   load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
  cellType=2;ub = 0.5; lb= 0.4; 
%   iidx = 1:length(data_nls);
%   mcellid = iidx%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
%  cids = data(cellType).cellIDs;
  
if(cellType==2)
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/predictions_cond3_4_OFF.mat');
else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/predictions_cond3_4_ON.mat');
end

  cellID =3318;%49
  iidx=1:length(cellID_log);
  idx=iidx(cellID_log==cellID);
 
% load true responses .

cols='rkrkrkrkrkrkkrkrkrkr';

       
condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
    igain_cnt=0;

for gain=gain_list
    igain_cnt=igain_cnt+1;
        R2_pl=[];
    icnt=1;pred=cell(2,1);clear ss;nTrials=30;
    for icond=3:4
        pred{icnt} = pred_log{idx,igain_cnt,icnt};    
        ss(icnt)=spkCondColl{icond};
        
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        convolve=15;
        binSz = 1/120;len=1200;
        [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        pred{icnt}(pred{icnt}>100000)=100000;
        [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});

        R2_pl =squeeze(R2_log(idx,igain_cnt,:));
        PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
        
       
        p{1} =pred{icnt};
  
        icnt=icnt+1;
        
    end
    
h=figure('units','normalized','Position',[1 1 0.35 0.25]);
PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
title(sprintf('gain: %d, Correlations: WN: %0.04f, Null:%0.04f',gain,R2_pl(1),R2_pl(2)));
 %xlim([2.5,10]);axis off
end  



