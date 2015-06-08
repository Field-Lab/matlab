% make figure. 
cellID = 3736;
folder = '/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/cell3736';

load('/Volumes/Analysis/2011-12-13-2/subunits/data008_2/conepreprocess.mat');
coneidx=1:size(datarun.cones.weights,1);
iidx=1:length(datarun.cell_ids);
ref_cell=iidx(datarun.cell_ids == cellID)%1:10;% 4 is interesting.
thr=0.3;          
%% 
            
spks=datarun.spike_rate(ref_cell,:);
conewts = (datarun.cones.weights(:,ref_cell) );
conewts_logical = conewts >= max(conewts)*thr;
cones =  coneidx(conewts_logical);
   nc=length(cones);

   % sta
        X =double(datarun.cone_inputs(:,conewts_logical));

nc=length(cones);

% calculate STA; 
idx=1:size(X,1);
spk_tm = idx(double(spks)>0);
sta=zeros(nc,6); Y=X';
for ispktm=spk_tm
    if(ispktm>6)
    sta=sta+Y(:,ispktm-5:ispktm)*double(spks(ispktm));
    end
end
sta=sta/double(sum(spks));
%%

figure('Color','w');

for isuu =1:5
subplot(2,5,isuu+5);
nFilters = nc - 4 +isuu-1; 
load(strcat(folder,sprintf('/quadGMLM_%d_su_%d.mat',cellID,nFilters)));
 
% make su_log
su_log=zeros(nc,nc);
for ifit=1:50
    fitGMLM = fitGMLM_log{ifit};
k_est = zeros(nFilters,nc);
    for isu=1:nFilters
    k_est(isu,:)=fitGMLM.Linear.filter{isu}(1:nc);
    end
  
    su=double((-k_est)==repmat(max((-k_est),[],1),[nFilters,1]));
su_log = su_log + su'*su;
end



    icnt=0;
    for isubunits=1:size(su,1)
        if(sum(su(isubunits,:))~=0)
            icnt=icnt+1;
            cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
        end
    end
    
    
    
  % make figure
       col='rgbcmykmrgb'
     
    
       
 
       for icone = 1:nc
            for jcone=1:nc
                iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                if(su_log(icone,jcone)>0 & icone~=jcone)
            plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone));
            hold on;
                end
            end
       end
       
           for icone = 1:nc
            for jcone=1:nc
                iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                if(su_log(icone,jcone)>0 & icone~=jcone)
            x1 = mean(datarun.cones.centers(pair,1));
            y1 = mean(-datarun.cones.centers(pair,2));
            text(x1,y1,sprintf('%d',su_log(icone,jcone)));
            hold on;
                end
            end
           end
       
       for icone=1:nc
        scatter(datarun.cones.centers(cones(icone),1) ,-datarun.cones.centers(cones(icone),2),abs(sta(icone,end-1))*2000,'filled','r');
         hold on;
       end
     xlim([min(datarun.cones.centers(cones,1))-2,max(datarun.cones.centers(cones,1))+2]);
     
     ylim([min(-datarun.cones.centers(cones,2))-2,max(-datarun.cones.centers(cones,2))+2])
       title(sprintf('n SU: %d',nFilters));
end


%% 
subplot(2,5,[2])
load(strcat(folder,sprintf('/quadGMLM_%d_interaction_alex.mat',cellID)));
sig_log=plot_cone_interaction(binnedResponsesbigd,maskedMov,paircones,false);
title(sprintf('data: %d w.r.t. %d',paircones(1),paircones(2)));

subplot(2,5,3);
hist(sig_log);
title(sprintf('sigma hist: %d',paircones(1)));



subplot(2,5,[4])
load(strcat(folder,sprintf('/quadGMLM_%d_interaction_alex.mat',cellID)));
paircones = paircones(end:-1:1);
sig_log=plot_cone_interaction(binnedResponsesbigd,maskedMov,paircones,false);
title(sprintf('data: %d w.r.t. %d',paircones(1),paircones(2)));

subplot(2,5,5);
hist(sig_log);
title(sprintf('sigma hist: %d',paircones(1)));
