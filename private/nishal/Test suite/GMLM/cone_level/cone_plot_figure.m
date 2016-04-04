
folder = '/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm/gmlm_detailed';
load('/Volumes/Analysis/2011-12-13-2/subunits/data008/conepreprocess.mat');
coneidx=1:size(datarun.cones.weights,1);
thr=0.33;
  %%
    
% make figure. 
for cellID =  [6332,4681,2851,5851,4460,3346,1157,2341,4306,4081];
folder = '/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/quad_gmlm/gmlm_detailed';
load('/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36/conepreprocess.mat');
coneidx=1:size(datarun.cones.weights,1);
thr=0.5; %[827,1157,1907,2341,2851,3346,3961,4081,4306,4460,4637,4681,4862,4876,4936,4998,5011,5176,5191,5251,5476,5701,5851,5911,6076,6151,6241,6332,7036];
%        
iidx=1:length(datarun.cell_ids);
ref_cell=iidx(datarun.cell_ids == cellID)%1:10;% 4 is interesting.
         

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

% h=figure('Color','w','PaperSize',[35,5],'PaperPosition',[0 0 35 5]);
% 
% for isuu =1:5
% subplot(1,5,isuu);
% nFilters = nc - 4 +isuu-1; 
% load(strcat(folder,sprintf('/quad_gmlm_%d_su_%d.mat',cellID,nFilters)));
% %load(strcat(folder,sprintf('/nnmf_%d_su_%d.mat',cellID,nFilters)));
%   
% % make su_log
% su_log=zeros(nc,nc);
% for ifit=1:50
%     fitGMLM = fitGMLM_log{ifit};
% k_est = zeros(nFilters,nc);
%     for isu=1:nFilters
%     k_est(isu,:)=fitGMLM.Linear.filter{isu}(1:nc);
%     end
%   
%     su=double((-k_est)==repmat(max((-k_est),[],1),[nFilters,1]));
% su_log = su_log + su'*su;
% end
% 
%     icnt=0;
%     for isubunits=1:size(su,1)
%         if(sum(su(isubunits,:))~=0)
%             icnt=icnt+1;
%             cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
%         end
%     end
%     
%     
%     
%   % make figure
%        col='rgbcmykmrgb'
%      
%     
%        
%  
%        for icone = 1:nc
%             for jcone=1:nc
%                 iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
%                 if(su_log(icone,jcone)>0 & icone~=jcone)
%             plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone)/2);
%             hold on;
%                 end
%             end
%        end
%        
%            for icone = 1:nc
%             for jcone=1:nc
%                 iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
%                 if(su_log(icone,jcone)>0 & icone~=jcone)
%             x1 = mean(datarun.cones.centers(pair,1));
%             y1 = mean(-datarun.cones.centers(pair,2));
%             text(x1,y1,sprintf('%d',su_log(icone,jcone)));
%             hold on;
%                 end
%             end
%            end
%        
%        for icone=1:nc
%         scatter(datarun.cones.centers(cones(icone),1) ,-datarun.cones.centers(cones(icone),2),abs(sta(icone,end-1))*2000/2,'filled','r');
%          hold on;
%        end
%      xlim([min(datarun.cones.centers(cones,1))-2,max(datarun.cones.centers(cones,1))+2]);
%      
%      ylim([min(-datarun.cones.centers(cones,2))-2,max(-datarun.cones.centers(cones,2))+2])
%        title(sprintf('n SU: %d',nFilters));
% end
%  print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/quad_gmlm/gmlm_summary/gmlm_%d.pdf',cellID))

%% TOP sub-units
h=figure('Color','w','PaperSize',[35,5],'PaperPosition',[0 0 35 5]);

for isuu =1:5
subplot(1,7,isuu);
nFilters = nc - 4 +isuu-1; 
if(nFilters<1)
    continue;
end
load(strcat(folder,sprintf('/quad_gmlm_%d_su_%d.mat',cellID,nFilters)));
  

% make SU log
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

       nsu_use = cell_su_list(cell_su_list(:,1)==cellID,2);
       su_wt_sort = sort(su_log(~logical(eye(nc,nc))),'descend');
       thr = su_wt_sort(nsu_use*2);
       
       
       
    icnt=0;
    for isubunits=1:size(su,1)
        if(sum(su(isubunits,:))~=0)
            icnt=icnt+1;
            cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
        end
    end
    
    
    
  % make figure
       %col='rgbcmykmrgb'
        col='g';
        su_col='b';
       for icone = 1:nc
            for jcone=1:nc
                iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                if(su_log(icone,jcone)>0 & icone~=jcone)
                    if(su_log(icone,jcone)>=thr)
           plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone)/2,'Color',su_col);
            hold on;
                    else
            %plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone)/2,'Color',col);
            hold on;    
                    end
                    
                end
            end
       end
       
  
       
%        for icone = 1:nc
%             for jcone=1:nc
%                 iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
%                 if(su_log(icone,jcone)>0 & icone~=jcone)
%             x1 = mean(datarun.cones.centers(pair,1));
%             y1 = mean(-datarun.cones.centers(pair,2));
%             text(x1,y1,sprintf('%d',su_log(icone,jcone)));
%             hold on;
%                 end
%             end
%         end
       
       for icone=1:nc
        scatter(datarun.cones.centers(cones(icone),1) ,-datarun.cones.centers(cones(icone),2),abs(sta(icone,end-1))*2000/2,'filled','r');
         hold on;
         if(nFilters==nc)
         text(datarun.cones.centers(cones(icone),1) ,-datarun.cones.centers(cones(icone),2),sprintf('%d',icone));
         end
       end
      
     xlim([min(datarun.cones.centers(cones,1))-2,max(datarun.cones.centers(cones,1))+2]);
     
     ylim([min(-datarun.cones.centers(cones,2))-2,max(-datarun.cones.centers(cones,2))+2])
       title(sprintf('nSU: %d, select %d',nFilters,nsu_use));
end

% s=hgexport('readstyle','cone_1x7');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/nnmf/nnmf_summary/nnmf_%d.eps',cellID),s);
 print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/quad_gmlm/gmlm_top/gmlm_%d.pdf',cellID))



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