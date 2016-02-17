
%% Initial stuff
%%
% Condition strings
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
location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

%% Load Movies
%  


% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies_OFF_addivitiy
    %[stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-2/Visual/%d.rawMovie',imov),rawMovFrames,1);
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_from_01/%d.rawMovie',imov),rawMovFrames,1);
    
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
for itime=3%[1,100]
subplot(3,1,1);
    imagesc(condMov{3}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot(3,1,2);
    imagesc(condMov{4}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot (3,1,3);
    imagesc(abs(condMov{3}(:,:,itime)' - condMov{4}(:,:,itime)'));axis image;colormap gray;colorbar;caxis([0,1]);
pause(1/120);
end


%% Predict response  and compare with real response! 
% 
% 
% % set up stuff for reading real responses
% dataRuns = dataRuns_OFF_additivity;
% 
% WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002';
% 
% datarun=load_data(WN_datafile)
% datarun=load_params(datarun)
% condDuration=10;
% nConditions=1;
% cellTypeId = 1;
% InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
% 
% for cellID =2328% InterestingCell_vis_id;
% cols='rkrkrkrkrkrkkrkrkrkr';
% spkCondColl=cell(length(dataRuns),1);
% close all;
% 
%     for idata=1:length(dataRuns)
%     Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data0%02d',dataRuns(idata));
%     neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
%     [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
%     end
% 
% % predict response
% %%
% 
% % destination = 'pc2015_11_09_1_analysis_fits/SUs_data025';
% % cellData001 = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)]);
% icnt=0;
% 
% 
% pred1=cell(2,1);PSTH_rec_log=cell(2,1);PSTH_pred_log=cell(2,1);
% clear('ss');
%  R2_pl = [];
% 
% for nSU=[1,2,3,4,5];        
% destination = 'pc2015_11_09_1_analysis_fits/SUs_data025_quad';
% cellData = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)])
% fitGMLM = cellData.fitGMLM_log{nSU};
%     
% scale_tf=1;
% 
% for icond=[2]%[1,2,3,4,5];
% icnt=icnt+1;
%     realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/1200,12000);
%         convolve=150;
%     binSz = 1/1200;len=12000;
% [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
% 
% 
% movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
% movd_masked = movd(repmat(cellData.mask>0,[1,1,size(movd,3)]));
% movd=movd;%*0.4804/sqrt(var(movd_masked(:)));
%  maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
%  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
% 
% nTrials=100;
% % [pred1{icnt},lam]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,1);
%   [pred1{icnt},lam]= predictGMLM_gamma2(fitGMLM,maskedMov2,nTrials,2,1);
%  pred1{icnt} = pred1{icnt}';
%  ss(icnt)=spkCondColl{icond};
%  
%   predd1{1} = pred1{icnt};
%   
%  [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,predd1{1});
%  %PSTH_pred=lam';
%  R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
%  PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
% end
% %close all;
% end
% 
% PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
% hh=plot_record_prediction(ss,pred1,'PSTH_struct',PSTH_str);
% title(sprintf('%0.02f,%0.02f,%0.02f,%0.02f,%0.02f',R2_pl(1),R2_pl(2),R2_pl(3),R2_pl(4),R2_pl(5)));
% %xlim([8,10]);
% %      if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s',datarun.cell_types{cellTypeId}.name)))
% %          mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s',datarun.cell_types{cellTypeId}.name));
% %      end
% %     s=hgexport('readstyle','raster');
% %     hgexport(hh,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,cellID),s);
% end

    %%
    %% prediction using sub-unit model
    load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
  cellType=2;
  ub = 2%0.5; 
  lb= 0%0.4; 
 movies = movies_OFF_addivitiy;


% load nls
 if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end


  iidx = 1:length(data_nls);
  mcellid = iidx;%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 cids = data(cellType).cellIDs(mcellid);
  
  
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
gain_list = [1,1.5,2];
pred_log = cell(length(mcellid),5,length(gain_list),2);
R2_log = zeros(length(mcellid),5,length(gain_list),2);
cellID_log=[];

 for imcell=1:length(mcellid)
     display(sprintf('Finished: %d out of %d cells',imcell,length(mcellid)));
  current_cell=mcellid(imcell);
  cellID = cids(imcell);
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
    


       cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2/Cell_%d',cellID));
   for    nSU=1:5;
       fitGMLM = cellData.fitGMLM_log{nSU};
    close all
    
    igain_cnt=0;

for gain=gain_list
    igain_cnt=igain_cnt+1;
        R2_pl=[];
    icnt=1;pred=cell(2,1);clear ss;nTrials=30;
    for icond=3:4
        
         movd = gain*(condMov{icond}-0.5);  % *2 if it is a low contrast movie!
         maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
         %  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
        [pred{icnt},lam]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,1);
        pred{icnt}=pred{icnt}';
%        pred{icnt} = %response_lnp(data_nls(current_cell),gain*(condMov{icond}-0.5),nTrials);    
        ss(icnt)=spkCondColl{icond};
        
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        convolve=15;
        binSz = 1/120;len=1200;
        [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        pred{icnt}(pred{icnt}>100000)=100000;
        [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});

        R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
        PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
        
        R2_log(imcell,nSU,igain_cnt,icnt)=R2_pl(end);
        p{1} =pred{icnt};
        pred_log(imcell,nSU,igain_cnt,icnt) =p;
        icnt=icnt+1;
        
    end
    
h=figure;%('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
title(sprintf('gain: %d, %0.04f,%0.04f',gain,R2_pl(1),R2_pl(2)));
% xlim([8,10]);
% pause
end  

   end
 end
%save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF.mat','pred_log','R2_log','gain_list','cellID_log');
 
%% 

load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');


cellType=2;
ub = 0.9%0.5; 
lb= 0.4;



iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
mcellid = iidx(cell_bin);
cids = data(cellType).cellIDs(cell_bin);

if(cellType==2)
     load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF_lam_max80.mat');
    %load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF_no_op_nl.mat');
else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_ON_no_op_nl.mat');
end
gain_list = [1,1.5,2];

figure;
icnt=0;
for nSU=1:5
for iplot=1:3
    icnt=icnt+1;
subplot(5,3,(icnt));
plot(R2_log(cell_bin,nSU,iplot,1),R2_log(cell_bin,nSU,iplot,2),'.');
xlim([0,1]);ylim([0,1]);
hold on;plot([0,1],[0,1],'g');
axis square
title(sprintf('nSU: %d, gain: %0.2f',nSU,gain_list(iplot)));
end
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
    %load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF_lam_max80.mat');
     %load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF_no_op_nl.mat');
     load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU_metrics.mat')

else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_ON.mat');
end

  cellID = 243%49%3318;
  iidx=1:length(cellID_log);
  idx=iidx(cellID_log==cellID);
 
% load true responses .

cols='rkrkrkrkrkrkkrkrkrkr';
gain_list = 1;%[1,1.5,2];
       
condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
for nSU=1:5
    igain_cnt=0;
for gain=gain_list
    igain_cnt=igain_cnt+1;
        R2_pl=[];
    icnt=1;pred=cell(2,1);clear ss;nTrials=30;
    for icond=3:4
        pred{icnt} = pred_log{idx,nSU,icnt};    
        ss(icnt)=spkCondColl{icond};
        
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        convolve=7;
        binSz = 1/120;len=1200;
        [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
       % pred{icnt}(pred{icnt}>100000)=100000;
        [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});

        R2_pl =squeeze(R2_log(idx,nSU,:));
        PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
        
       
        p{1} =pred{icnt};
  
        icnt=icnt+1;
        
    end
    
h=figure('units','normalized','Position',[1 1 0.35 0.25]);
PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
hh=plot_record_prediction(ss,pred);%,'PSTH_struct',PSTH_str);
title(sprintf('#SU: %d, gain: %0.02f,Correlation : WN :%0.04f, Null: %0.04f',nSU,gain,R2_pl(1),R2_pl(2)));
% xlim([2.5,10]);axis off
end  

end

%%
    %% prediction using output non-linearities
    load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
   figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
  cellType=2;
  ub = 0.5; 
  lb= 0.4; 
 movies = movies_OFF_addivitiy;


% load nls
 if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end


  iidx = 1:length(data_nls);
  mcellid = iidx;%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 cids = data(cellType).cellIDs(mcellid);
  
  
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
gain_list = [1,1.5,2];
pred_log = cell(length(mcellid),5,length(gain_list),2);
R2_log = zeros(length(mcellid),5,length(gain_list),2);
cellID_log=[];

 for imcell=1:length(mcellid)
     display(sprintf('Finished: %d out of %d cells',imcell,length(mcellid)));
  current_cell=mcellid(imcell);
  cellID = cids(imcell);
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
    


       cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2/Cell_%d',cellID));
       cellData_opnl =  load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl/Cell_%d',cellID));
       
   for    nSU=1:5;
       fitGMLM = cellData.fitGMLM_log{nSU};
       fitGMLM.NL = cellData_opnl.fit_nl(nSU);
    close all
    
    igain_cnt=0;
    
    for gain=gain_list
        igain_cnt=igain_cnt+1;
        R2_pl=[];
        icnt=1;pred=cell(2,1);clear ss;nTrials=30;
        for icond=3:4
            
            movd = gain*(condMov{icond}-0.5);  % *2 if it is a low contrast movie!
            maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
            %  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
            
            [pred{icnt},lam]= predictGMLM_bias_lr_opnl(fitGMLM,maskedMov,nTrials,1);
            pred{icnt}=pred{icnt}';
            %        pred{icnt} = %response_lnp(data_nls(current_cell),gain*(condMov{icond}-0.5),nTrials);
            ss(icnt)=spkCondColl{icond};
            
            realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
            convolve=15;
            binSz = 1/120;len=1200;
            [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
            %pred{icnt}(pred{icnt}>100000)=100000;
            [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});
            
            R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
            PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
            
            R2_log(imcell,nSU,igain_cnt,icnt)=R2_pl(end);
            p{1} =pred{icnt};
            pred_log(imcell,nSU,igain_cnt,icnt) =p;
            icnt=icnt+1;
            
        end
        
        h=figure;%('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
        PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
        hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
        title(sprintf('gain: %d, %0.04f,%0.04f',gain,R2_pl(1),R2_pl(2)));
        % xlim([8,10]);
        % pause
    end
    
   end
 end
save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_OFF_lam_max80.mat','pred_log','R2_log','gain_list','cellID_log');
 
    %% prediction using optimal output non-linearity to match the response
    
    load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
    figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
    dataRuns = dataRuns_OFF_additivity;
    cellType=2;
    ub = 0.5;
    lb= 0.4;
    movies = movies_OFF_addivitiy;


% load nls
 if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end


  iidx = 1:length(data_nls);
  mcellid = iidx;%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 cids = data(cellType).cellIDs(mcellid);
  
  
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
cond_list=[3,4,5,6];
% iterate over different cells

sus_list = [1,2,3,4,5,6,7,8,9,10];
nSUs= length(sus_list);
pred_log = cell(length(mcellid),nSUs,length(cond_list));
R2_log = zeros(length(mcellid),nSUs,length(cond_list));
cellID_log=[];
nl_info_log=cell(length(mcellid),nSUs,length(cond_list));
error_cells=[];
metrics  = zeros(length(mcellid),nSUs,length(cond_list),4);

for imcell=1:length(mcellid)
    pause
    display(sprintf('Starting: %d out of %d cells',imcell,length(mcellid)));
    current_cell=mcellid(imcell);
    cellID = cids(imcell);
    % load true responses .
    cellID_log=[cellID_log;cellID];
    if(data_nls(current_cell).cellID== cellID)
        display('cell IDs match');
    end
    cols='rkrkrkrkrkrkkrkrkrkr';
    try 
    condDuration=10;
    nConditions=1;
    for idata=1:length(dataRuns)
        Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
        neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
        [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
    if cellType==2
    cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2/Cell_%d',cellID));
    else
    cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/ON Parasol2/Cell_%d',cellID));
  
    end
    
    for    nSU=sus_list%1:nSUs
        nSU
        fitGMLM = cellData.fitGMLM_log{nSU};
        close all
        R2_pl=[];
        icnt=1;pred=cell(2,1);clear ss;nTrials=30;
        for icond=[3,4,5,6]
            % real response
            realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
            convolve=7;
            binSz = 1/120;len=1200;
            [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
            ss(icnt)=spkCondColl{icond};
            [nSU,icond,sum(cellData.mask(:))]
            % predicted response
            [real_resp_bin1,time]=calculate_psth_fcn2(1,binSz,len,realResp);
            movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
            maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
            [pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_bias_lr_optimal_op_nl(fitGMLM,maskedMov,nTrials,1,real_resp_bin1);
            %[pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_gamma2_lr_optimal_op_nl(fitGMLM,maskedMov,nTrials,1,real_resp_bin1);
            %[pred{icnt},lam]= predictGMLM_gamma2_lr(fitGMLM,maskedMov,nTrials,2,1); pred{icnt}=pred{icnt}';
            %[pred{icnt},lam]= predictGMLM_full(fitGMLM,maskedMov,nTrials); pred{icnt}=pred{icnt}';
         
            nl_info_log{imcell,nSU,icnt} = 0;
            
            % [pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,1);
            % pred{icnt} = pred{icnt}';
            % pause
           
            ss(icnt)=spkCondColl{icond};
            [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});
           
            R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
            PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
            
            R2_log(imcell,nSU,icnt)=R2_pl(end);
            
            metrics(imcell,nSU,icnt,:) = compute_raster_metrics(realResp,pred{icnt});
            p{1} =pred{icnt};
            pred_log(imcell,nSU,icnt) =p;
            icnt=icnt+1;
            
        end
        
        h=figure;%('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
        PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
        hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
        title(sprintf(' #SU:%d , %0.04f,%0.04f,%0.04f,%0.04f',nSU,R2_pl(1),R2_pl(2),R2_pl(3),R2_pl(4)));
        %title(sprintf(' %0.04f,%0.04f,%0.04f,%0.04f',R2_pl(1)));
        % xlim([8,10]);
         pause(1)
    end
     catch
         display(sprintf('Error in cell: %d, conditions: %d',imcell,icond));
         error_cells = [error_cells,cids(imcell)];
         
     end
    
end
%save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_OFF_metrics_optimalNL.mat','metrics','pred_log','R2_log','error_cells','cellID_log','cond_list','nl_info_log','sus_list');

%% 
icond=3;icnt=1;
realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
nTrials = 6;%6;
m = compute_raster_metrics (realResp(1:nTrials,:),pred{icnt}(1:nTrials,:));
[m(1),m(3)]

%% see results
load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
cellType=2;%ub = 0.6; lb= 0.4;
lb=0.4;ub=0.6;
if(cellType==1)
    data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_ON_optimalNL_nSU10.mat')
else
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_OFF_optimalNL_nSU10.mat')
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU.mat')
    data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_OFF_metrics_optimalNL.mat')
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU_metrics.mat')
end

iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
mcellid = iidx(cell_bin);
cids = data(cellType).cellIDs(cell_bin);

figure;
icnt=0;
for isu=data_pred.sus_list
    icnt=icnt+1;
subplot(2,5,icnt)
plot(data_pred.R2_log(cell_bin,isu,1),data_pred.R2_log(cell_bin,isu,2),'.','markersize',5);
hold on;
plot([0,1],[0,1],'g')
xlim([0,1]);ylim([0,1]);
axis square
title(sprintf('nSU: %d ',isu));
end

figure;
for isu=data_pred.sus_list

    subplot(2,5,isu)
    for icond =[1,2,3,4]
    plot(data_pred.metrics(cell_bin,isu,icond,3),data_pred.metrics(cell_bin,isu,icond,1),'.');  
    hold on;
    end
    %legend('WN','Null','WN+Null var controlled','WN+Null var not controlled');
    title(sprintf('# SU %d:',isu));
end

% vary trials
clear vary_trials
icond=3;icnt=1;isu=5;

m_log=[]; 
for itrial=1:5
vary_trials(itrial).m_log=[];
end

for icell = mcellid
    cellID = data(cellType).cellIDs(icell);
    cellID
    condDuration=10;
    nConditions=1;
    close all
    for idata=1:length(dataRuns)
        Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
        neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
        [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    itrial=0;
    for nTrials = [4,6,8,12,20]
        itrial=itrial+1;
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        pred=data_pred.pred_log(icell,isu,icnt);
        if(~isempty(pred{1}))
        m = compute_raster_metrics(realResp(1:floor(30/nTrials):nTrials*floor(30/nTrials),:),pred{1}(1:floor(30/nTrials):nTrials*floor(30/nTrials),:));
        %m_log=[m_log;[m(3),m(1)]];
        else
            m=[0,0,0,0];
        end
        vary_trials(itrial).m_log=[vary_trials(itrial).m_log;[m(3),m(1)]];
    end
    
end

figure;
for itrial = 1:length(vary_trials)
    hold on;
    plot(vary_trials(itrial).m_log(:,1),vary_trials(itrial).m_log(:,2),'k.','markersize',15);
end

for icnt =[1,2]
    plot(data_pred.metrics(cell_bin,isu,icnt,3),data_pred.metrics(cell_bin,isu,icnt,1),'.','markersize',15,'Color',cols(icnt,:));
    hold on;
end


%legend('WN','Null','WN+Null var controlled','WN+Null var not controlled');
title(sprintf('# SU %d:',isu));



figure;
for icond=1:2;
cols = distinguishable_colors(2);
for iisu=1:length(data_pred.sus_list(1:end-1))
    isu = data_pred.sus_list(iisu);
    jsu = data_pred.sus_list(iisu+1);
    metric_ori = data_pred.metrics(cell_bin,isu,icond,1)./data_pred.metrics(cell_bin,isu,icond,3);
    metric_end = data_pred.metrics(cell_bin,jsu,icond,1)./data_pred.metrics(cell_bin,jsu,icond,3);
    for icell=1:sum(cell_bin)
        plot([isu+(icond-1)*0.1,jsu+(icond-1)*0.1],[metric_ori(icell),metric_end(icell)],'-*','Color',cols(icond,:));
        hold on;
    end
end
end

