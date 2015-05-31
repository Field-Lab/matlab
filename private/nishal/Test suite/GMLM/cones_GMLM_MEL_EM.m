       %% Load alex cone data
            load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');
           close all
           % load('/Volumes/Analysis/2011-12-13-2/subunits/data008_2/conepreprocess.mat');
           % % Good cell - WN , cell ID: 166 thr=0.3
            
           %  load('/Volumes/Analysis/2011-12-13-2/subunits/data009/conepreprocess.mat');
           % Good cell - Voronoi, cell ID : 168 ,thr =0.5
            
           

            coneidx=1:size(datarun.cones.weights,1);

            ref_cell=2;% 4 is interesting.
            thr=0.3;
            cellID=datarun.cell_ids(ref_cell)
           
            
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
            nFilters =nc-1 % on average, 2 sub-units!
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
%%    
    binnedResponsesbigd = double(spks(1,2:end)');
    maskedMov = X(1:end-1,:)';       
  %%
%             binnedResponsesbigd = double(spks(1,1:end)');
%             tf=sta(maxCone,end:-1:1);
%             maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
%             % maskedMov = maskedMov(:,2:end); spks = spks(1:end-1);
            %%
            sd = sqrt(diag(maskedMov*maskedMov'/size(maskedMov,2)));
            maskedMov = maskedMov.*repmat(1./sd,[1,size(maskedMov,2)]);
            
            figure;hist(maskedMov(:));title('Filtered Movie')
            maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

            
            interval=10;
         
 
 mov_use=maskedMov;
 filteredStimDim=size(mov_use,1);
 nSU = nFilters;
 
 k_est_log=cell(50,1);
 

su_log=zeros(nc,nc);
fitGMLM_log=cell(50,1);
 for ifit=1:50
 % [fitGMLM,output] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
%
% 
% spike.home 
%  [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,spike.home,mov_use);
%  
% figure;
%  plot(fitGMLM_full2.hist.hexpanded)
 
 
 %
    k_est = zeros(nFilters,nc);
    for isu=1:nFilters
    k_est(isu,:)=fitGMLM.Linear.filter{isu}(1:nc);
    end
    
    k_est_log{ifit}=k_est;

k_est=k_est_log{ifit};
    su=double((-k_est)==repmat(max((-k_est)),[nFilters,1]));
su_log = su_log + su'*su;
fitGMLM_log{ifit} = fitGMLM;
end
%su_log=su_log/50;
save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/Cell%d.mat',cellID) , 'fitGMLM_log');

    icnt=0;
    for isubunits=1:size(su,1)
        if(sum(su(isubunits,:))~=0)
            icnt=icnt+1;
            cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
        end
    end


       
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
       
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/quadGMLM_%d_su_1.eps',cellID));
%close all

       %% SU activation plot 
       su=[2,3];
       su_activation_plot_gamma(fitGMLM,mov_use,2,su);
       iso_response_gamma(binnedResponsesbigd,mov_use,fitGMLM,2,su)
       %% Predict
       R2_log=zeros(50,1);
       for ifit=1:50
           ifit
       fitGMLM=fitGMLM_log{ifit};
       gamma=2;nTrials=1; 
       rec=binnedResponsesbigd; mov_pred=mov_use;
       [pred_trials,lam] = predictGMLM_gamma2(fitGMLM,mov_pred,nTrials,gamma,interval); % interva1=10;
      %[pred,lam] = predictGMLM_bias(fitGMLM,mov_pred,nTrials,interval);
    
        pred=lam; 
     
      pred_ss = zeros(length(pred)/10,1);
       for itime=1:length(pred_ss)
       pred_ss(itime) = sum(pred((itime-1)*10+1:(itime)*10));
       end
       
       xxx=[pred_ss,rec];
       norm(double(rec) -double(pred_ss))/norm(double(pred_ss));
      
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
       end