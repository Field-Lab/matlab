function  [WNSTA,h,h2,dimensions]= spike_triggred_analysis_pc2015_03_09_3(condMov,mov_bin_sz,neuronPath,cellID,nConditions,condDuration,condNoMov,WN_datafile_short,d_save,GLM_fit_link,Mask,dimensions_in)
%% Take movie
movie = condMov{condNoMov}-0.5;
movieLen=size(movie,3);
%% Load spikes and triggers
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=neuronFile.getSpikeTimes(cellID);
TTL=double(neuronFile.getTTLTimes());

spks=double(CellSpkTimes);

rawMovFrames = condDuration*120;
TTLperCondperTrial=floor(rawMovFrames/100)+1;

nTrials=floor(length(TTL)/(TTLperCondperTrial*nConditions));

spkColl=cell(nTrials,1);

%condDuration=12;
% nConditions=6;
samplesPerCondition=condDuration*20000;
ConditionStartTimes=[0:nConditions]*samplesPerCondition+1;

spkCondColl=struct('spksColl',[]);
spkCondColl(nConditions).spksColl=zeros(nTrials,movieLen);
TTLstartOld=0;

for iTrial=1:nTrials

TTLstart=TTL((iTrial-1)*TTLperCondperTrial*nConditions+1);
TTLend=TTLstart+condDuration*nConditions*20000;%TTL((iTrial)*TTLperTrial+1);
(TTLstart-TTLstartOld)/20000
spkColl{iTrial}=(spks((spks>=TTLstart)&(spks<=TTLend))'-TTLstart);


for icond=1:2
    
TTLstart=TTL(((iTrial-1)*nConditions+icond-1)*TTLperCondperTrial+1);
TTLend=TTLstart+condDuration*20000;%TTL((iTrial)*TTLperTrial+1);
(TTLstart-TTLstartOld)/20000

% 
% spkCondColl(icond).spksColl{iTrial}=(spks((spks>=TTLstart)&(spks<=TTLend))'-TTLstart);
% if(isempty(spkCondColl(icond).spksColl{iTrial}))
% spkCondColl(icond).spksColl{iTrial}=1;
% end

for iframes = 1:movieLen
tstart =TTL(((iTrial-1)*nConditions+icond-1)*TTLperCondperTrial+1 +  floor(iframes/100)) + rem(iframes-1,100)*mov_bin_sz*20000;
tend = tstart+mov_bin_sz*20000; 
spkCondColl(icond).spksColl(iTrial,iframes)=sum(((spks>=tstart)&(spks<=tend)))~=0;
end
TTLstartOld=TTLstart;
end
end

% Have spikes for the condition

%% Calculate STA
  movie4D = zeros(size(movie,1),size(movie,2),3,movieLen);
    for iframe=1:size(movie,3);
    movie4D(:,:,1,iframe)=movie(:,:,iframe);
    movie4D(:,:,2,iframe)=movie(:,:,iframe);
    movie4D(:,:,3,iframe)=movie(:,:,iframe);
    end
    mov_params.mov=movie4D;
    
    sta_params.Filtlen=30;

    cell_params.binsPerFrame=1;
    
    response.spksGen=spkCondColl(condNoMov).spksColl;
    aa=repmat([1:movieLen],[1,1]);
    response.mov_frame_number=aa(:);
    
    WNSTA_coll=cell(nTrials,1);
    WNSTA=zeros(size(movie,1),size(movie,2),sta_params.Filtlen);
    for iTrial=1:nTrials
    sta_params.useTrial=iTrial; % Need to use all 5 trials
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
    WNSTA_coll{iTrial} = response.analyse.STA;
    WNSTA = WNSTA+WNSTA_coll{iTrial};
    end
    WNSTA=WNSTA/nTrials;
    
         figure
         for itime=1:sta_params.Filtlen
             itime
         imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
         caxis([min(WNSTA(:)),max(WNSTA(:))]);
         colorbar
         pause(1/120)
         end
         
         % Clip STA
         cell_p.STAlen=sta_params.Filtlen;
         stas{1}=zeros(size(WNSTA,1),size(WNSTA,2),3,cell_p.STAlen);
         for itime = 1:cell_p.STAlen
         stas{1}(:,:,1,itime)=WNSTA(:,:,itime);
         stas{1}(:,:,2,itime)=WNSTA(:,:,itime);
         stas{1}(:,:,3,itime)=WNSTA(:,:,itime);
         end
         
         [new_stas,totalMaskAccept,CellMasks]=clipSTAs(stas,cell_p);
            if(sum(Mask(:))==0)
                Mask = CellMasks{1};
            end
       
       
       STA_flat=zeros(size(WNSTA,1)*size(WNSTA,2),sta_params.Filtlen);
          for itime=1:sta_params.Filtlen
             itime
         xx=((WNSTA(:,:,itime)).*Mask)';
         STA_flat(:,itime)=xx(:);
          end
          
          h=figure('Color','w');
         subplot(5,4,[13,14,15,16,17,18,19,20]);
          [U,S,V]=svds(STA_flat,20);
         plot(diag(S),'*');
         title('Re-STA spectrum');

        subplot(5,4,[5,6,7,8,9,10,11,12]);
        fit_glm = load(strcat(GLM_fit_link,sprintf('/%d.mat',cellID)));
        
         average_time=mean(STA_flat,1);
         average_time=average_time -mean(average_time);
         plot((average_time)/max(abs(average_time)),'b');
         hold on;
         x= fit_glm.fittedGLM.linearfilters.Stimulus.time_rk1;
         x=x-mean(x);
         plot(-(x)/max(abs(x)),'r');
         l=legend('Null STA - avg','Orignial GLM','Location','best');
         set(l,'FontSize',8);
         ylim([-1,1]);
         title('Re-STA Temporal filter');
            
         
           [v,i]=max(abs(average_time));
           subplot(5,4,[1,2]);
           imagesc(WNSTA(:,:,i)');
           colormap gray
           caxis([min(WNSTA(:)),max(WNSTA(:))]);
           colorbar
           axis image
           title('Strongest STA frame');
           
          subplot(5,4,[3,4]);
          imagesc(Mask');
          
%% FIT GLM , fixed linear front end / rank 1?
% movie_full=zeros(size(movie,1),size(movie,2),nTrials*size(movie,3));
% spikes_rec=zeros(size(spkCondColl(condNoMov).spksColl,2),1);
% 
% for itrial=1:nTrials
% spikes_rec((itrial-1)*movieLen+1:itrial*movieLen)=spkCondColl(condNoMov).spksColl(itrial,:);
% movie_full(:,:,(itrial-1)*movieLen+1:itrial*movieLen)=movie;
% end
% idx=1:length(spikes_rec);
% spk_tms=idx(spikes_rec==1);
% movie_xml = 'BW-8-4-0.48-11111';
% 
% try
% fittedGLM=glm_fit_from_movie_spikes({cellID},WN_datafile_short, movie_xml, nTrials*condDuration, d_save,movie_full,spk_tms,WNSTA);
% catch
% load(strcat(d_save,sprintf('/%d.mat',cellID))); 
% end

%% STC? 

reSTA = WNSTA;

Filtdim1=size(reSTA,1);
Filtdim2=size(reSTA,2);
Filtlen=size(reSTA,3);
reSTA = reSTA(logical(repmat(Mask,[1,1,Filtlen])));

%xxreSTA = reshape(reSTA,[Filtdim1*Filtdim2*Filtlen,1]);
xxreSTA=reSTA;
xxreSTA=xxreSTA/norm(xxreSTA(:));

reSTC=zeros(length(xxreSTA));
indx=[59:1:movieLen];
binnedResponses = spkCondColl(condNoMov).spksColl;
for iTrial=1:nTrials

binnedResponsesTrial=binnedResponses(iTrial,:);

framesValid = indx(binnedResponsesTrial(indx)>0);

for iframe=framesValid
    iframe
    a = movie(:,:,iframe:-1:iframe-Filtlen+1);
    a=a(logical(repmat(Mask,[1,1,Filtlen])));
  xx= a;
%  mask!

%  xx=  reshape(mov_new2(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
% Don't mask!

xxnew = (xx-((xx'*xxreSTA)*xxreSTA));
  
reSTC=reSTC+ xxnew*xxnew'*binnedResponsesTrial(iframe);
end
end
reSTC=reSTC / (sum(binnedResponses(:))-1);

[u,s,v]=svds(reSTC,100);
figure;
plot(diag(s),'*')

%% plot points 
h2=figure('Color','w');
try
    xxreSTA=dimensions_in.xxreSTA;
    u=dimensions_in.u;
    Mask=dimensions_in.Mask;
    s=dimensions_in.s;
    display('Values Loaded')
catch

end
component1 = xxreSTA;
component2 = u(:,1);
subplot(3,2,1);
proj_plot_script_pc2015_03_09_3
title('STA - STC 1');

component1 = xxreSTA;
component2 = u(:,2);
subplot(3,2,2);
proj_plot_script_pc2015_03_09_3
title('STA - STC 2');

component1 = u(:,1);
component2 = u(:,2);
subplot(3,2,3);
proj_plot_script_pc2015_03_09_3
title('STC 1 - STC 2');

component1 = u(:,3);
component2 = u(:,2);
subplot(3,2,4);
proj_plot_script_pc2015_03_09_3
title('STC 3 - STC 2');

component1 = u(:,1);
component2 = u(:,4);
subplot(3,2,5);
proj_plot_script_pc2015_03_09_3
title('STC 1 - STC 4');

subplot(3,2,6);
plot(diag(s),'*');
title('STC magnitudes');

%%
dimensions.xxreSTA=xxreSTA;
dimensions.u=u;
dimensions.Mask=Mask;
dimensions.s=s;
end