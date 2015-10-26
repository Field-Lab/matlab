%% prediction to null stimulus for GLM fits. 
% Dataset: 2015-03-09-2/data038 (WN) and data042 (null run), NDF0

%% Load dataset information

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
wrong_xml = 'BW-8-2-0.48-22222-40x40';
stim_length=1800;% 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);



%% Load different null movies
nConditions=6;
% make movies
interval=2;
condMov=cell(nConditions,1);
rawMovFrames=1272/(2);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
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
            condMov{icnt}(:,:,ifcnt)=double(qq(:,:,iframe)); % cond mov is between 0 and 1 now!
        end
        
    end
    
end

%% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];



%% Do predictions - full

cell_list = datarun.cell_types{2}.cell_ids;

 cellID=1531

data1 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full.mat',cellID,cellID));
data2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
user_STA_depth=30;
extract_movie_response2
%  if(~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID),'dir'))
%         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID));
%  end
% print(hhh,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/cellID_%d_mask.pdf',cellID,cellID));

fitGMLM_log=cell(5,1);
for isu=1:3
fitGMLM_log{isu} = data1.fitGMLM_full2_log{isu};
end
for isu=4:5
fitGMLM_log{isu} = data2.fitGMLM_full2_log{isu};
end
mask=data2.totalMaskAccept2;
%%
nSU=5
fitGMLM=fitGMLM_log{nSU};
   % fitGMLM = scramble_fitGMLM_mix(fitGMLM);

   
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
%% multiplicative gains 

gain_range = [1:0.2:6,6.5,7,7.5,8,8.5,9];

    convolve=150;
    len = 12720;
    binSz=1/1200;
icond=5;%1:nConditions
spkMat = makeSpikeMat(spkCondColl(icond).spksColl,binSz,len);
[PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,spkMat);
jgainrange = [0.1,0.3,0.6,0.8,1,2,2.2,2.4,2.6,2.8,3];

fr = zeros(length(jgainrange),length(gain_range)); var_log=zeros(length(jgainrange),length(gain_range)); R2_log = var_log ; 
fr_l=fr; var_l = var_log; R2_l = R2_log;
nfits =1 ;
nl = @(x)x;
for ifit =1:nfits
    
fr = zeros(length(jgainrange),length(gain_range)); var_log=zeros(length(jgainrange),length(gain_range)); R2_log = var_log ; 
jidx=0;
for jgain = jgainrange
    jidx = jidx+1
    iidx=0;
for igain=gain_range
    iidx=iidx+1;
close all;
% Make predictions
pred1=cell(1,1);

movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=29;
 [pred1{1},lam]= predictGMLM_full_gains(fitGMLM,maskedMov,nTrials,nl,jgain);
 pred1{1} = pred1{1}';
 fr(jidx,iidx) =sum(pred1{1}(:))/(size(pred1{1},1) * size(pred1{1},2)/1200);
  
 % Calculate psth
 [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred1{1});
    
    % Recorded 
 var_log(jidx,iidx) = var(PSTH_pred);
 
   R2_log(jidx,iidx) = R_2_value(PSTH_rec',PSTH_pred'); 

end
end
R2_l = R2_l + R2_log;
var_l = var_l + var_log;
fr_l = fr_l + fr;
end

R2_l=R2_l/nfits;
var_l = var_l/nfits;
fr_l = fr_l/nfits;

fr_rec = spkCondColl(icond).avgSpkRate
var_rec = var(PSTH_rec)
figure;
subplot(1,4,1);
[C_fr,h] = contourf(gain_range,jgainrange,fr_l,[fr_rec,fr_rec]);
subplot(1,4,2);
[C_var,h]=contourf(gain_range,jgainrange,var_l,[var_rec,var_rec]);
subplot(1,4,3);
contourf(gain_range,jgainrange,R2_l,20);

hh=figure;
subplot(1,4,1);
[~,h] = contourf(gain_range,jgainrange,fr_l,20);

subplot(1,4,2);
[~,h]=contourf(gain_range,jgainrange,var_l,20);

subplot(1,4,3);
contourf(gain_range,jgainrange,R2_l,20);

subplot(1,4,4)
plot(C_fr(1,2:end),C_fr(2,2:end));hold on;
plot(C_var(1,2:end),C_var(2,2:end));

%legend('Observed average firing rate','Observed PSTH varaince');
%%

figure;
plot(gain_range,fr);
hold on;
plot(gain_range,spkCondColl(icond).avgSpkRate *ones(size(gain_range)));

%%
realResp = makeSpikeMat(spkCondColl(icond).spksColl,1/120,1272);
fr= mean(realResp,1)/(1/120);

igain=1.9%2%2%input('What gain to choose?');
jgain = 0.7%0.8%0.6275%0.6275;
pred1=cell(1,1);
icond=5%5%1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 R2_pl = [];
 for ipred = 1:20
nTrials=100;
 [pred1{1},lam]= predictGMLM_full_gains(fitGMLM,maskedMov,nTrials,nl,jgain);
 pred1{1} = pred1{1}';
 [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred1{1});
 R2_pl(ipred) = R_2_value(PSTH_rec',PSTH_pred')
 end
 R2_p = mean(R2_pl)
 
   h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
    % close all
    %figure;
    ss = spkCondColl(icond);
   hh=plot_record_prediction(ss,pred1);

   %% multiplicative for stimulus, and power law of o/p NL
jgainrange = [0.005,0.01,0.05,0.1,0.3,0.6,0.8,1,2,2.2,2.4,2.6,2.8,3];
gain_range =  [1:0.2:6,6.5,7,7.5,8,8.5,9];

    convolve=150;
    len = 12720;
    binSz=1/1200;
icond=5;%1:nConditions
spkMat = makeSpikeMat(spkCondColl(icond).spksColl,binSz,len);
[PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,spkMat);


fr = zeros(length(jgainrange),length(gain_range)); var_log=zeros(length(jgainrange),length(gain_range)); R2_log = var_log ; 
fr_l=fr; var_l = var_log; R2_l = R2_log;
nfits =1 ;

for ifit =1:nfits
    
fr = zeros(length(jgainrange),length(gain_range)); var_log=zeros(length(jgainrange),length(gain_range)); R2_log = var_log ; 
jidx=0;
for jgain = jgainrange
    nl = @(x)(x.^jgain);
    jidx = jidx+1
    iidx=0;
for igain=gain_range
    iidx=iidx+1;
close all;
% Make predictions
pred1=cell(1,1);

movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=29;
 [pred1{1},lam]= predictGMLM_full_gains(fitGMLM,maskedMov,nTrials,nl,1);
 pred1{1} = pred1{1}';
 fr(jidx,iidx) =sum(pred1{1}(:))/(size(pred1{1},1) * size(pred1{1},2)/1200);
  
 % Calculate psth
 [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred1{1});
    
    % Recorded 
 var_log(jidx,iidx) = var(PSTH_pred);
 
   R2_log(jidx,iidx) = R_2_value(PSTH_rec',PSTH_pred'); 

end
end
R2_l = R2_l + R2_log;
var_l = var_l + var_log;
fr_l = fr_l + fr;
end

R2_l=R2_l/nfits;
var_l = var_l/nfits;
fr_l = fr_l/nfits;

fr_rec = spkCondColl(icond).avgSpkRate
var_rec = var(PSTH_rec)
figure;
subplot(1,4,1);
[C_fr,h] = contourf(gain_range,jgainrange,fr_l,[fr_rec,fr_rec]);
subplot(1,4,2);
[C_var,h]=contourf(gain_range,jgainrange,var_l,[var_rec,var_rec]);
subplot(1,4,3);
contourf(gain_range,jgainrange,R2_l,20);

hh=figure;
subplot(1,4,1);
[~,h] = contourf(gain_range,jgainrange,fr_l,20);

subplot(1,4,2);
[~,h]=contourf(gain_range,jgainrange,var_l,20);

subplot(1,4,3);
contourf(gain_range,jgainrange,R2_l,20);

subplot(1,4,4)
plot(C_fr(1,2:end),C_fr(2,2:end));hold on;
plot(C_var(1,2:end),C_var(2,2:end));

%legend('Observed average firing rate','Observed PSTH varaince');
%%

realResp = makeSpikeMat(spkCondColl(icond).spksColl,1/120,1272);
fr= mean(realResp,1)/(1/120);

igain=3%2%2%input('What gain to choose?');
jgain =0.7%0.8%0.6275%0.6275;
nl = @(x) (x.^jgain);
pred1=cell(1,1);
icond=5%5%1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 R2_pl = [];
 for ipred = 1:1
nTrials=100;
 [pred1{1},lam]= predictGMLM_full_gains(fitGMLM,maskedMov,nTrials,nl,1);
 pred1{1} = pred1{1}';
 [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred1{1});
 R2_pl(ipred) = R_2_value(PSTH_rec',PSTH_pred')
 end
 R2_p = mean(R2_pl)
 
   h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
    % close all
    %figure;
    ss = spkCondColl(icond);
   hh=plot_record_prediction(ss,pred1);

    hh=plot_record_pred12(ss,pred4,pred1);
    
