addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
%% Dataset details
WN_datafile = '2015-03-09-2/d18-28-norefit/data026-from-d18-28/data026-from-d18-28';
WN_datafile_short='2015-03-09-2/d18-28-norefit/data026-from-d18-28/data026-from-d18-28';
movie_xml = 'BW-4-2-0.48-11111-80x80';
stim_length=1800;% in seconds

%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

for cellID =3721%datarun.cell_types{2}.cell_ids;
    cellID
%cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%%
extract_movie_response2;

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen)*0.8)];
testData=[floor(length(spksGen)*0.8)+1 :length(spksGen)];
 binnedResponsesbigd = spksGen(trainData);
 mov_use=maskedMovdd(:,trainData);
 
 filteredStimDim=size(mov_use,1); 

 for nSU =1:15;
 binnedResponsesbigd = spksGen(trainData);
 mov_use=maskedMovdd(:,trainData);
 
 for ifit=1:1
% [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 fitGMLM_sulog{nSU}=fitGMLM;
 %[fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
% [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,spike.home,mov_use);
% figure;
% plot(fitGMLM_full2.hist.hexpanded)
 end
 % Training Err 
  nTrials=1;
  binnedResponsesbigd_train = spksGen(trainData);
 mov_use_train=maskedMovdd(:,trainData);
 [pred_train{nSU},lam_train{nSU}]= predictGMLM_bias(fitGMLM,mov_use_train,nTrials,interval); 
 rec=binnedResponsesbigd_train;
 [f_val_train(nSU),R2_train(nSU)]=f_r2_test(lam_train{nSU},rec,interval);
 
 % Testing.
 nTrials=1;
  binnedResponsesbigd_test = spksGen(testData);
 mov_use_test=maskedMovdd(:,testData);
 [pred_test{nSU},lam_test{nSU}]= predictGMLM_bias(fitGMLM,mov_use_test,nTrials,interval); 
 rec=binnedResponsesbigd_test;
 [f_val_test(nSU),R2_test(nSU)]=f_r2_test(lam_test{nSU},rec,interval);

 end
% save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/Cell%d_full',cellID),'fitGMLM_log','fitGMLM_full2','mov_use','binnedResponsesbigd','nSU','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
 
figure('Color','w');
plot(1:length(f_val_train),f_val_train,'--*');
hold on;
plot(1:length(f_val_test),f_val_test,'--*');
legend('Training','Testing');
xlabel('number of SU');
ylabel('negative Likelihood');
title('Training and Testing Likelihood value');

figure('Color','w');
plot(1:length(R2_train),R2_train,'--*');
hold on;
plot(1:length(R2_test),R2_test,'--*');
legend('Training','Testing');
xlabel('number of SU');
ylabel('R2');
title('Training and Testing R-sq value');
end
 
%% Prediction 
 
  %% Show learned filters;
  mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(dim1,dim2,nSU);

figure;
for ifilt=1:nSU
subplot(2,2,ifilt)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.1);
xx=xx(x_coord,y_coord);
figure;
imagesc(xx);

%iso_response_bias_gmlm(binnedResponsesbigd,mov_use,fitGMLM);

%su_activation_plot(fitGMLM_full2,mov_use);
%% Response prediction 

%% Load different movies - NSEM
nMovies=30;

movie_total = zeros(320,160,nMovies*120);
for ichunk=1:nMovies
    ichunk
    a=load(sprintf( '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_%d.mat',ichunk));
    movie_total(:,:,(ichunk-1)*120 + 1 :ichunk*120)=a.movie;
end
movie_total = movie_total/255-0.5; % orientation unclear! 
movie_total=movie_total(81:240,:,:);
movie_len = size(movie_total,3);

movie_ld = zeros(dim1,dim2,movie_len);
sub_sample = 2; % stixel size of STA or model / 2
for ix = 1:dim1
    for iy=1:dim2
        movie_ld(ix,iy,:)=squeeze(mean(mean(movie_total((ix-1)*sub_sample +1 : ix*sub_sample,(iy-1)*sub_sample +1 : iy*sub_sample,:),1),2));
    end
end

%% WN repeats


%% Load recorded data
% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/d18-28-norefit/data024-from-d18-28';
neuronPath = [Null_datafile,sprintf('/data024-from-d18-28.neurons')];
condDuration=30;
nConditions=1;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);

%% Do predictions - full
% Make predictions
predf=cell(nConditions,1);

movd = movie_ld;
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=20;
 predf{1}= predictGMLM_full(fitGMLM_full2,maskedMov,nTrials)';


plot_record_prediction(spkCondColl,predf)

%% Do predictions - no feedback

% Make predictions
pred=cell(nConditions,1); lam=cell(nConditions,1);
for  icond=1:nConditions
movd = movie_ld;
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=20;
interval=1;
 [pred{icond},lam{icond}]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,interval);
 pred{icond}=pred{icond}';
  %pred{icond}= predictGMLM(fitGMLM,maskedMov,nTrials)';
 % pred{icond}= predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,2)';
end

plot_record_prediction(spkCondColl,pred)      

plot_record_prediction3(spkCondColl,pred1cell,pred10cell);


%% R-sq value
for icond=1:6
       predd=lam{icond}; 
     
       pred_ss = zeros(length(predd)/10,1);
       for itime=1:length(pred_ss)
       pred_ss(itime) = sum(predd((itime-1)*10+1:(itime)*10));
       end
       
       pred_ss=repmat(pred_ss,[29,1]);
       
       rec = makeSpikeMat(spkCondColl(icond).spksColl,1/120,condDuration/(1/120));
       rec=rec';
       rec=rec(:);
       
       
       % R2 value method 2
       x1 = pred_ss; y1 = rec; n=length(y1);
       r = (n*x1'*y1 - sum(x1)*sum(y1))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y1.^2) - sum(y1)^2));
       R2_log(icond) = r^2;
       
end
R2_log

%% Load WN movie
movie_xml = 'BW-4-2-0.48-11111-80x80';
stim_length=29.9;
stim_description=movie_xml;
xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
dashes=find(stim_description=='-');
StimulusPars.type=stim_description(1:dashes(1)-1);
StimulusPars.pixelsize = str2double(stim_description(dashes(1)+1:dashes(2)-1));
StimulusPars.refreshrate = str2double(stim_description(dashes(2)+1:dashes(3)-1));
StimulusPars.RNG = str2double(stim_description(dashes(3)+1:dashes(4)-1));
% try
%     StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
%     StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
% catch
%     StimulusPars.height = 32; 
%     StimulusPars.width = 64;
% end
StimulusPars.tstim = 1/120;
fitframes=stim_length*120; % seconds * 120 frames per second / interval


disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, [0:100/120:stim_length], fitframes/StimulusPars.refreshrate);
temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
fitmovie_color=zeros(width,height,3,fitframes);

try
    StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
    StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
catch
    StimulusPars.height =height; 
    StimulusPars.width = width;
end


for i=1:fitframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/StimulusPars.refreshrate));
end

clear temp_fitmovie i 

movie_ld=squeeze(mean(fitmovie_color,3))-0.5;

%% 
%% Load recorded data
% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/d18-28-norefit/data025-from-d18-28';
neuronPath = [Null_datafile,sprintf('/data025-from-d18-28.neurons')];
condDuration=29.9;
nConditions=1;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


%%

lam1=lam_test{1};
lam2=lam_test{5};

tdiff = abs(lam1-lam2)>0.02;
tdiff2=repmat(tdiff,[10,1]);
tdiff_fill=zeros(length(tdiff),1);
for itime =1:10:length(tdiff_fill)
tdiff_fill(itime:itime+9) = sum(tdiff(itime:itime+9))>0;
end
tdiff_fill=logical(tdiff_fill);


[f_val_test1,R2_test1]=f_r2_test(lam1(tdiff_fill),rec(tdiff_fill(1:10:end)),interval)
[f_val_test2,R2_test2]=f_r2_test(lam2(tdiff_fill),rec(tdiff_fill(1:10:end)),interval)