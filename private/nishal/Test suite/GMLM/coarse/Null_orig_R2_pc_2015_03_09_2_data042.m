%% Do predictions - full

cell_list =1531;

for cellID=cell_list

R2_log=cell(6,1);
PSTH_rec=cell(6,1);
for icond=1:6
R2_log{icond} = zeros(17,50);

    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
  spkMat = makeSpikeMat(spkCondColl(icond).spksColl,1/120,1272);
    [PSTH_rec{icond},time]=calculate_psth_fcn2(50,1/120,1272,spkMat);
end
for nSU=1:17
    nSU
data = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/OFF parasol/detailed/CellID_%d/exp_gmlm/Cell%d_quad_gmlm_su_%d.mat',cellID,cellID,nSU));
mask=data.totalMaskAccept;

    
    %fitGMLM = scramble_fitGMLM_mix(fitGMLM);
igain=1;
    for ifit=1:50
        ifit
        close all
        
        fitGMLM=data.fitGMLM_log{ifit};

% Make predictions
pred1=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(data.ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
interval=1;
 [pred1{icond},lam]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,interval);

     PSTH_pred = lam;
    
    % Recorded 
    R2_log{icond}(nSU,ifit) = R_2_value(PSTH_rec{icond}',lam);
end

    end
end

end

%% Calculate PSTH and find correlation between recorded and predictions.

icond=1;
r2_su = R2_log{icond};
icond=5;
r2_su_null = R2_log{icond};
nc=17;
       h=figure;
      
       errorbar(1:nc,mean(r2_su,2),sqrt(var(r2_su')),'-*');
       hold on;
       plot(1:nc,max(r2_su'),'-*');
       hold on;
       errorbar(1:nc,mean(r2_su_null,2),sqrt(var(r2_su_null')),'-*');
       hold on;
       plot(1:nc,max(r2_su_null'),'-*');
       hold on;
     
       
       legend('mean for WN','max for WN','mean for null','max for null','location','best');
       title('R^2 value v/s max. number of sub-units');
       xlabel('Number of sub-units');
       ylabel('R^2');
       hgexport(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/')