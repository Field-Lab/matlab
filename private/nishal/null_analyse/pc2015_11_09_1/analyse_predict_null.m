
%% Initial stuff
%%
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
%% Load Movies
%  


% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies_ON_additivity
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
for itime=3%[1,100]
subplot(3,1,1);
    imagesc(condMov{3}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot(3,1,2);
    imagesc(condMov{4}(:,:,itime)');axis image;colormap gray;caxis([0,1]);colorbar
subplot (3,1,3);
    imagesc(abs(condMov{3}(:,:,itime)' - condMov{4}(:,:,itime)'));axis image;colormap gray;colorbar;caxis([0,1]);
pause(1/120);
end
%% show frames 


% set up stuff for reading real responses
dataRuns = dataRuns_ON_additivity;
WN_datafile ='/Volumes/Analysis/2015-11-09-1/streamed/data009/data009'

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
condDuration=10;
nConditions=1;
cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 


% Display stimulus
figure;
for itime=150%[1,100]
for icond=1:6
    subplot(3,2,icond);
    imagesc(repelem(condMov{icond}(:,end:-1:1,itime),20,20)'-0.5);axis image;colormap gray;caxis([-0.5,+0.5]);colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
    
  plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
 end
end

% Display stimulus - magnified
figure;
for itime=150%[1,100]
for icond=1:6
    subplot(3,2,icond);
    imagesc(repelem(condMov{icond}(:,end:-1:1,itime),20,20)'-0.5);axis image;colormap gray;caxis([-0.5,0.5]);colorbar;set(gca,'xTick',[]);set(gca,'yTick',[]);
    
  plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
xlim([500,700]);
ylim([200,400]);
end
end

% get mask 
% mask = (cMap{2}-cMap{1})~=0;
% imagesc(mask);

figure; 
for icond=1:6
subplot(3,2,icond);
hist(condMov{icond}(repelem(mask,1,1,3))-0.5,15);
[x,N] = hist(condMov{icond}(repelem(mask,1,1,3))-0.5,15);
hold on
title(sprintf('cond: %d',icond));
xlim([-0.5,0.5]);
set(gca,'yTick',[]);

if(icond==1 || icond==2)
    hold on;
plot([0.4804,0.4804],[0,1.1*max(x)],'r'); hold on;
plot([-0.4804,-0.4804],[0,1.1*max(x)],'r'); 
end

if(icond>=3)
    hold on;
plot([0.2412,0.2412],[0,1.1*max(x)],'r'); hold on;
plot([-0.2412,-0.2412],[0,1.1*max(x)],'r'); 

end
set(gca,'xTick',[-0.48,-0.24,0.24,0.48]);
end

%% Predict response  and compare with real response! 


% set up stuff for reading real responses
dataRuns = dataRuns_ON_additivity;
WN_datafile ='/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data009/data009'

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
condDuration=10;
nConditions=1;
cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 

for cellID =2328% InterestingCell_vis_id;
cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
close all;

    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-11-09-1/d14_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end

% predict response
%%

% destination = 'pc2015_11_09_1_analysis_fits/SUs_data025';
% cellData001 = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)]);
icnt=0;


pred1=cell(2,1);PSTH_rec_log=cell(2,1);PSTH_pred_log=cell(2,1);
clear('ss');
 R2_pl = [];

for nSU=[1,2,3,4,5];        
destination = 'pc2015_11_09_1_analysis_fits/SUs_data025_quad';
cellData = load(['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)])
fitGMLM = cellData.fitGMLM_log{nSU};
    
scale_tf=1;

for icond=[2]%[1,2,3,4,5];
icnt=icnt+1;
    realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/1200,12000);
        convolve=150;
    binSz = 1/1200;len=12000;
[PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);


movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
movd_masked = movd(repmat(cellData.mask>0,[1,1,size(movd,3)]));
movd=movd;%*0.4804/sqrt(var(movd_masked(:)));
 maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

nTrials=100;
% [pred1{icnt},lam]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,1);
  [pred1{icnt},lam]= predictGMLM_gamma2(fitGMLM,maskedMov2,nTrials,2,1);
 pred1{icnt} = pred1{icnt}';
 ss(icnt)=spkCondColl{icond};
 
  predd1{1} = pred1{icnt};
  
 [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,predd1{1});
 %PSTH_pred=lam';
 R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
 PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
end
%close all;
end

PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
hh=plot_record_prediction(ss,pred1,'PSTH_struct',PSTH_str);
title(sprintf('%0.02f,%0.02f,%0.02f,%0.02f,%0.02f',R2_pl(1),R2_pl(2),R2_pl(3),R2_pl(4),R2_pl(5)));
%xlim([8,10]);
%      if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s',datarun.cell_types{cellTypeId}.name)))
%          mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s',datarun.cell_types{cellTypeId}.name));
%      end
%     s=hgexport('readstyle','raster');
%     hgexport(hh,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_11_09_1/d_pred_null_using_longnull_100tr/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,cellID),s);
    end
