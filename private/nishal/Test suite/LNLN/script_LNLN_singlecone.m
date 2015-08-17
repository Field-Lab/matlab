%% Single cone level script

%% Make model 
model = model_LNLN_midget();

%% Generate test response
% sz = model.gridSzX/ (1*3);
% 
% Tlen = 120*10;
% movie = gpuArray((randn(sz,sz,Tlen)>0)-0.5);
% dt=10/120;
% nTrials=30;
% [response,~] = generateResp_LNLN_singlecone(model,movie*1.5,dt,nTrials);
% h= figure;
% plotSpikeRaster(response~=0,'PlotType','vertline');
% %print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/sample_firing.pdf'));
% 
% Tlen = 120*60*3*10/10;
% movie = gpuArray((randn(sz,sz,Tlen)>0)-0.5);
% dt=10/120;
% nTrials=1;
% [response,~] = generateResp_LNLN_singlecone(model,movie*5,dt,nTrials);

% %% Calulate STA 
% idx = 1:Tlen;
% spktm = idx(response~=0 & idx >4);
% 
% STA= zeros(sz,sz,3);
% 
% for itime=spktm
% STA = STA + response(itime)*movie(:,:,itime-2:itime);
% end
% 
% STA = STA / sum(response);
% STA = gather(STA);
% 
% figure;
% for itime=1:3
%     itime
% imagesc(STA(:,:,itime));
% colormap gray
% caxis([min(STA(:)),max(STA(:))]);
% pause
% end
% 
% % STA strongest frame with cone map;
% strongestFrame = STA(:,:,26);
% szstr = size(strongestFrame,1);
% 
% ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);
% 
% h=figure;
% imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*30 + model.totalConeMap3D);
% axis image
% title('STA 30 min');
% set(gca,'xTick',[]);
% set(gca,'yTick',[]);
% print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA30.pdf'));
% 
% h=figure;
% imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*20 + model.totalConeMap3D);
% axis image
% title('STA 20 min');
% set(gca,'xTick',[]);
% set(gca,'yTick',[]);
% print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA20.pdf'));

%% Directly give cone inputs...


Tlen = 120*30;
coneinp_movie = gpuArray(randn(model.nCones,Tlen)*10);
dt=1/120;
nTrials=30;
[response,~] = generateResp_LNLN_cone_inps(model,coneinp_movie,dt,nTrials);

h= figure;
plotSpikeRaster(response~=0,'PlotType','vertline');


Tlen = 120*60*30;
coneinp_movie = gpuArray(randn(model.nCones,Tlen)*10);
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN_cone_inps(model,coneinp_movie,dt,nTrials);


% calculate STA - just one frame! 
idx = 1:Tlen;
spktm = idx(response~=0);

STA= zeros(model.nCones,1);

for itime=spktm
STA = STA + response(itime)*coneinp_movie(:,itime);
end

STA = STA / sum(response);
STA = gather(STA);

% plot STA 
STA2D = zeros(model.gridSzX,model.gridSzY);
for icone=1:model.nCones
STA2D = STA2D + model.cone_data{icone}.coneMap*STA(icone);
end

STA2D = STA2D(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);
figure;
imagesc(STA2D);
colormap gray
set(gca,'xTick',[]);
set(gca,'yTick',[]);


%% Run model.
filteredStimDim = size(coneinp_movie,1);
for nFrontEnds = 8%[9:-1:1,10];
nFrontEnds
interval=1;
gamma=2;
su_log =zeros(model.nCones);
fitGMLM_log =cell(50,1);
f_val_log = [];
for ifit=1:50
ifit
[fitGLM,f_val] = fitGMLM_MEL_EM_power2(response',gather(coneinp_movie),filteredStimDim,nFrontEnds,interval,gamma);

fitGMLM_log{ifit} = fitGLM;
f_val_log(ifit) = f_val;
% show model fits

k_est = zeros(nFrontEnds,model.nCones);
for isu=1:nFrontEnds
    k_est(isu,:)=fitGLM.Linear.filter{isu}(1:model.nCones);
end
su=double((k_est)==repmat(max((k_est),[],1),[nFrontEnds,1]));
su_log = su_log + su'*su

end
save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d.mat',nFrontEnds),'fitGMLM_log','f_val_log','su_log');

end

%%
cols = distinguishable_colors(nFrontEnds+1);
figure;
uall = zeros(model.gridSzX,model.gridSzY,3);

for isu=1:nFrontEnds

 u =zeros(model.gridSzX,model.gridSzY);
for icone = 1:model.nCones
u =u+model.cone_data{icone}.coneMap*su(isu,icone)*10;
end
u = repelem(u,1,1,3);
u(:,:,1) = cols(isu+1,1)*u(:,:,1);
u(:,:,2) = cols(isu+1,2)*u(:,:,2);
u(:,:,3) = cols(isu+1,3)*u(:,:,3);

uall = uall+u;


end
uall = uall(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);
close all
figure
imagesc(uall);
colormap gray;
pause

%% R2 value ? 
% test movie and response

Tlen = 120*30;
coneinp_movie_test = gpuArray(randn(model.nCones,Tlen)*10);
dt=1/120;
nTrials=30;
[resptest,~] = generateResp_LNLN_cone_inps(model,coneinp_movie_test,dt,nTrials);

R2_log=zeros(10,50);
for nSU =1:10
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d',nSU));
nSU
for ifit=1:50
fitGLM =fitGMLM_log{ifit};
interval = 1;
[predictedResponse,lam] = predictGMLM_gamma2_lr(fitGLM,gather(coneinp_movie_test),1,2,interval);
v = R_2_value(lam,resptest(1,:)');
R2_log(nSU,ifit) =v;

end
end

r2_su = R2_log;
nc=10;
       h=figure;
       %plot(1:nc,mean(r2_su_tr,2),'-*');
       hold on;
       errorbar(1:nc,mean(r2_su,2),sqrt(var(r2_su')),'-*');
       hold on;
       %plot(1:nc,max(r2_su_tr'),'-*');
       hold on;
       plot(1:nc,max(r2_su'),'-*');
       hold on;
      % plot(1:nc,min(r2_su_tr'),'-*');
      % hold on;
      % plot(1:nc,min(r2_su'),'-*');
       
       legend('mean across 50 fits','max across 50 fits','location','best');
       title('R^2 value v/s max. number of sub-units');
       xlabel('Number of sub-units');
       ylabel('R^2');
 hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_R2_su_test_individual_trail.eps'));

%% 
  
for nSU=1:10
 load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d',nSU));
    close all

    nc=model.nCones;
    h=figure;
            for icone = 1:nc
                for jcone=1:nc
                    iconeidx =icone; jconeidx=jcone; pair = [iconeidx;jconeidx];
                    if(su_log(icone,jcone)>0 & icone~=jcone)
                        plot(model.conesY(pair,1),-model.conesX(pair,1),'LineWidth',su_log(icone,jcone)/2);
                        hold on;
                    end
                end
            end
            
            for icone = 1:nc
                for jcone=1:nc
                    iconeidx =icone; jconeidx=jcone; pair = [iconeidx;jconeidx];
                    if(su_log(icone,jcone)>0 & icone~=jcone)
                        x1 = mean(model.conesY(pair,1));
                        y1 = mean(-model.conesX(pair,1));
                        text(x1,y1,sprintf('%d',su_log(icone,jcone)));
                        hold on;
                    end
                end
            end
            
            for icone=1  :nc
                scatter(model.conesY(icone,1) ,-model.conesX(icone,1),abs(STA(icone))*100,'filled','r');
                hold on;
            end
            xlim([min(model.conesY)-5,max(model.conesY)+5]);
            
            ylim([min(-model.conesX)-5,max(-model.conesX)+5])
            title(sprintf('max no. of SU: %d',nSU));
            
            set(gca,'xTick',[]);
            set(gca,'yTick',[]);
           hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d.eps',nSU));
end


%% 
nc=model.nCones;
  x=repmat(1:nc,[nc,1]);
  y=repmat([1:nc]',[1,nc]);
  
 mask = logical(x>y)
for nSU=1:7
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d.mat',nSU));
    close all
ss_l = su_log(mask);
h=figure;
hist(ss_l,25);
xlim([0,55]) 
set(gca,'yTick',[]);
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_hist_su_%d.eps',nSU));
pause
end