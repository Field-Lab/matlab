%% Load alex cone data
%  load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');
% close all

 load('/Volumes/Analysis/2011-12-13-2/subunits/data008_2/conepreprocess.mat');
  thr=0.3;
% Good cell - WN , cell ID: 166 thr=0.3

%load('/Volumes/Analysis/2011-12-13-2/subunits/data009/conepreprocess.mat');
%            Good cell - Voronoi, cell ID : 168 ,thr =0.5
%thr=0.5

% load('/Volumes/Analysis/2011-12-13-2/subunits/data008/conepreprocess.mat');
% thr = 0.33

%               load('/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36/conepreprocess.mat');
%               thr=0.5; %[827,1157,1907,2341,2851,3346,3961,4081,4306,4460,4637,4681,4862,4876,4936,4998,5011,5176,5191,5251,5476,5701,5851,5911,6076,6151,6241,6332,7036];
% %

coneidx=1:size(datarun.cones.weights,1);
iidx=1:length(datarun.cell_ids);
for cellID=3736%[6332,4681,2851,5851,4460,3346,1157,2341,4306,4081]%[827,1157,2341,2851,3346,3961,4081,4306,4460,4637,4681,4862,4876,4936,4998,5011,5176,5191,5251,5476,5701,5851,5911,6076,6151,6241,6332,7036];%1321%[1351,2251,3586,4576,5162,1321];
    ref_cell=iidx(datarun.cell_ids == cellID)%166%1:10;% 4 is interesting.
    
    %cellID=datarun.cell_ids(ref_cell);
    spks=datarun.spike_rate(ref_cell,:);
    
    %% Real cone wts
    conewts = (datarun.cones.weights(:,ref_cell) );
    conewts_logical = conewts >= max(conewts)*thr;
    cones =  coneidx(conewts_logical);
    %% Logical cone wts
    %             conewts_logical = logical(datarun.cones.weights(:,ref_cell) );
    %              cones = coneidx (conewts_logical);
    %%
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
    h=figure;
    plot(sta');
    
    [maxCone,time] = find(abs(sta)== max(max(abs(sta))));
    %
    %             title(sprintf('CellID: %d',cellID));
    %             hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/STA_%d.eps',cellID));
    %% for WN movie
    % 2nd from lat
        binnedResponsesbigd = double(spks(1,2:end)');
        maskedMov = X(1:end-1,:)';
    
        % 3rd from last
    %        binnedResponsesbigd = double(spks(1,3:end)');
    %     maskedMov = X(1:end-3,:)';
    %% for voronoi, filter the movie with STA , for NNMF
    %             binnedResponsesbigd = double(spks(1,1:end)');
    %             tf=mean(sta(:,end:-1:1),1);
    %             maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
    %             % maskedMov = maskedMov(:,2:end); spks = spks(1:end-1);
    
    %% for voronoi, filter the movie with STA , for GMLM-Quad
% %     binnedResponsesbigd = double(spks(1,1:end)');
% %     tf=-mean(sta(:,end:-1:1),1);
% %     maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
% %     % maskedMov = maskedMov(:,2:end); spks = spks(1:end-1);
    %%
    sd = sqrt(diag(maskedMov*maskedMov'/size(maskedMov,2)));
    maskedMov = maskedMov.*repmat(1./sd,[1,size(maskedMov,2)]);
    
    figure;hist(maskedMov(:));title('Filtered Movie')
    maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
    
    % Partition the data
    
    trainFrac = 0.8 %[0.2:-0.05:0.05,0.04:-0.002:0.01,0.008:-0.0002:0.002]%0.8:-0.2:0.2; % default by Jeremy
        subDivide = 10; % default by Jeremy
        
        
        %%%% Partition of the data
        [train test] = partitionDat(maskedMov',trainFrac,subDivide,1);
        train=train';test=test';
        
        %%%% Partition of the data
        [resptrain resptest] = partitionDat(binnedResponsesbigd,trainFrac,subDivide,1);
        
        
        
        %%
        
        interval=10;
        mov_use=train; %resptrain = binnedResponsesbigdtrain;
        filteredStimDim=size(mov_use,1);
        fitGMLM_su_log = cell(nc,1); 
        fval_su_log = cell(nc,1);
        su_log_su = cell(nc,1);
        
        for nFilters = [1:nc-5];%[nc-2:-1:nc-4,nc] %[nc-4: nc-0] % on average, 2 sub-units!
            nSU = nFilters;
            
            k_est_log=cell(50,1);
            
            
            su_log=zeros(nc,nc);
            fitGMLM_log=cell(50,1);
            f_val_log=[];
            for ifit=1:50
                ifit
                % [fitGMLM,output] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);
                 [fitGMLM,output] = fitGMLM_MEL_EM_power2(resptrain,mov_use,filteredStimDim,nSU,interval,2);
                %[fitGMLM,output] = fitGMLM_EM_power2(resptrain,mov_use,filteredStimDim,nSU,interval,2);
                %
                %
                % spike.home
                %  [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,spike.home,mov_use);
                %
                % figure;
                %  plot(fitGMLM_full2.hist.hexpanded)
                
                f_val_log=[f_val_log;output];
                %
                k_est = zeros(nFilters,nc);
                for isu=1:nFilters
                    k_est(isu,:)=fitGMLM.Linear.filter{isu}(1:nc);
                end
                
                k_est_log{ifit}=k_est;
                
                k_est=k_est_log{ifit};
                su=double((-k_est)==repmat(max((-k_est),[],1),[nFilters,1]));
                su_log = su_log + su'*su
                fitGMLM_log{ifit} = fitGMLM;
            end
            %su_log=su_log/50;
            %save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/Cell%d.mat',cellID) , 'fitGMLM_log');
            fitGMLM_su_log{nSU} = fitGMLM_log; 
            fval_su_log{nSU} = f_val_log;
            su_log_su{nSU} = su_log; 
            
            icnt=0;
            for isubunits=1:size(su,1)
                if(sum(su(isubunits,:))~=0)
                    icnt=icnt+1;
                    cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
                end
            end
            
            
            
            %
            col='rgbcmykmrgb'
            figure;
            for isu=1:length(cells(ref_cell).sub_units)
                for icone=cells(ref_cell).sub_units(isu).indices
                    plot(datarun.cones.centers(icone,1),-datarun.cones.centers(icone,2),strcat(col(isu),'*'));
                    axis equal
                    hold on
                end
                plot(datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,1) ,-datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,2),col(isu));
                hold on;
            end
            
            close all;
            
            h=figure;
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
            title(sprintf('Cell ID: %d',cellID));
            
            set(gca,'xTick',[]);
            set(gca,'yTick',[]);
            save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_%d_su_%d_data.mat',cellID,nFilters),'fitGMLM_log','f_val_log','su_log');
            hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_%d_su_%d.eps',cellID,nFilters));
        end
 
%save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_%d_data.mat',cellID),'fitGMLM_su_log','f_val_su_log','su_log_su');
                
end


%% Cone interaction 
paircones=[4,5]
% plot_cone_interaction(resptrain , mov_use, paircones);
plot_cone_interaction(binnedResponsesbigd,maskedMov,paircones);
save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/quadGMLM_%d_interaction_alex.mat',cellID),'binnedResponsesbigd','maskedMov','paircones');
       %% SU activation plot 
       su=[2,3];
       su_activation_plot_gamma(fitGMLM,mov_use,2,su);
       iso_response_gamma(resptrain ,mov_use,fitGMLM,2,su)
       %% Predict
       ll_su = zeros(nc,50);
       ll_su_tr = zeros(nc,50);
       
       for nSU = 1:nc
       loadedData =load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_3736_su_%d_data.mat',nSU));
       fitGMLM_log = loadedData.fitGMLM_log;
       f_val_log = loadedData.f_val_log;
       
       
       R2_log=zeros(50,1);
       ll=[];
       r2=[];
       for ifit=1:50
           ifit
       fitGMLM=fitGMLM_log{ifit};
       gamma=2;nTrials=1; interval=10;
       rec=resptest; mov_pred=test;
       [pred_trials,lam] = predictGMLM_gamma2_lr(fitGMLM,mov_pred,nTrials,gamma,interval); % interva1=10;
       if(nSU==1)
       lam1=lam;
       end
      %[pred,lam] = predictGMLM_bias(fitGMLM,mov_pred,nTrials,interval);
       ll(ifit) = (-sum(lam)*interval/120 + rec'*log(lam))/length(lam); 
       r2(ifit) = R_2_value(lam,rec);
       end
       ll_su(nSU,:) = -ll';
       r2_su(nSU,:) = r2';
       
       ll=[];
       r2=[];
       for ifit=1:50
           ifit
       fitGMLM=fitGMLM_log{ifit};
       gamma=2;nTrials=1; interval=10;
      rec=resptrain; mov_pred=train;
       [pred_trials,lam] = predictGMLM_gamma2_lr(fitGMLM,mov_pred,nTrials,gamma,interval); % interva1=10;
      %[pred,lam] = predictGMLM_bias(fitGMLM,mov_pred,nTrials,interval);
       ll(ifit) = (-sum(lam)*interval/120 + rec'*log(lam))/length(lam); 
       r2(ifit) = R_2_value(lam,rec);
       end
      
       ll_su_tr(nSU,:) = f_val_log';
       r2_su_tr(nSU,:) = r2';
       end
       
       
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
 hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data008_gmlm_fig/quad_gmlm_R2_3736_su_%d.eps',nSU));

      %%
         pred_ss = zeros(length(pred)/10,1);
       for itime=1:length(pred_ss)
       pred_ss(itime) = sum(pred((itime-1)*10+1:(itime)*10));
       end
       
       %xxx=[pred_ss,rec];
       %norm(double(rec) -double(pred_ss))/norm(double(pred_ss))
      
%        % R2 value - 
%        y=rec;
%        ybar = mean(y);
%        f=pred_ss;
%        SStot = sum((y-ybar).^2);
%        SSreg = sum((f-ybar).^2);
%        SSres = sum((y-f).^2);
%        R2 = 1 - (SSres/SStot)
%        
       
       % R2 value method 2
       x1 = pred_ss; y1 = rec; n=length(y1);
       r = (n*x1'*y1 - sum(x1)*sum(y1))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y1.^2) - sum(y1)^2));
       R2_log(ifit) = r^2;
%       end
       
       %% Moment based
       su_log2 = 1000*(train*train')/size(train,2)
       %
       col='rgbcmykm'
       figure;
       for isu=1:length(cells(ref_cell).sub_units)
           for icone=cells(ref_cell).sub_units(isu).indices
       plot(datarun.cones.centers(icone,1),-datarun.cones.centers(icone,2),strcat(col(isu),'*'));
       axis equal
       hold on
           end
           plot(datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,1) ,-datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,2),col(isu));
           hold on;
       end

       
       
      h=figure;
       for icone = 1:nc
            for jcone=1:nc
                iconeidx =cones(icone); jconeidx=cones(jcone); pair = [iconeidx;jconeidx];
                if(su_log2(icone,jcone)>0 & icone~=jcone)
            plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log2(icone,jcone));
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
       title(sprintf('Cell ID: %d',cellID));
       
       set(gca,'xTick',[]);
       set(gca,'yTick',[]);
       


%%
nc = model.nCones;
x=repmat(1:nc,[nc,1]);
y=repmat([1:nc]',[1,nc]);
  
 mask = logical(x>y)
for nSU=1:nc
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midget_su_%d.mat',nSU));
    close all
ss_l = su_log(mask);
h=figure;
hist(ss_l,25);
xlim([0,55]) 
set(gca,'yTick',[]);
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/model_cell/midgethist_su_%d.eps',nSU));
pause
end
       