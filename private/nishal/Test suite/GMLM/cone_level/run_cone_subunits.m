data = load('/Volumes/Analysis/2016-02-17-4/subunits/data_for_nishal_4876.mat');


nc = size(data.filt_inputs,1);
select_cones=[5,6,8] ;
su_log = zeros(numel(select_cones),numel(select_cones));
nSU=numel(select_cones);
interval=1;

for ifit=1:50
    ifit
    
    [fitGMLM,output] = fitGMLM_MEL_EM_power2(data.spike_rate,data.filt_inputs(select_cones,:),numel(select_cones),nSU,interval,2);
    
    k_est = zeros(nSU,numel(select_cones));
    for isu=1:nSU
        k_est(isu,:)=fitGMLM.Linear.filter{isu};
    end
    close all
    su=double((k_est)==repmat(max((k_est),[],1),[nSU,1]));
    su_log = su_log + su'*su;
    
end

%% plot estimated sub-units
figure;
for isu=1:nSU
    subplot(1,nSU,isu);
    for icone=1:nc
        
        scatter(centers((icone),1) ,-centers((icone),2),5,'filled','b');hold on;
        scatter(centers((icone),1) ,-centers((icone),2),max(k_est(icone,isu),0),'filled','r');
        axis equal
        hold on;
    end
    set(gca,'xTick',[]);set(gca,'yTick',[]);
xlim([min(centers(:,1))-2,max(centers(:,1))+2]); 
ylim([min(-centers(:,2))-2,max(-centers(:,2))+2])
end

%% plot map


centers = data.cone_coords;

h=figure;
for icone = 1:nc
    for jcone=1:nc
        pair = [icone;jcone];
        if(su_log(icone,jcone)>0 & icone~=jcone)
            plot(centers(pair,1),-centers(pair,2),'LineWidth',su_log(icone,jcone));
            hold on;
        end
    end
end

for icone = 1:nc
    for jcone=1:nc
        pair = [icone;jcone];
        if(su_log(icone,jcone)>0 & icone~=jcone)
            x1 = mean(centers(pair,1));
            y1 = mean(-centers(pair,2));
            text(x1,y1,sprintf('%d',su_log(icone,jcone)));
            hold on;
        end
    end
end

for icone=1:nc
    scatter(centers((icone),1) ,-centers((icone),2),20,'filled','r');
    hold on;
end
xlim([min(centers(:,1))-2,max(centers(:,1))+2]); 
ylim([min(-centers(:,2))-2,max(-centers(:,2))+2])
