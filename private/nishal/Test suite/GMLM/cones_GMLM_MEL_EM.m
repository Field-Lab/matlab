       %% Load alex cone data
            load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');
           close all
            %load('/Volumes/Analysis/2012-04-13-1/subunits/data002/conepreprocess.mat');

            coneidx=1:size(datarun.cones.weights,1);

for  ref_cell=1:10;% 4 is interesting.
            
            cellID=datarun.cell_ids(ref_cell)
           
            
            spks=datarun.spike_rate(ref_cell,:);
            X =double(datarun.cone_inputs(:,logical(datarun.cones.weights(:,ref_cell))));
            cones = coneidx (logical(datarun.cones.weights(:,ref_cell)));
            nc=length(cones);
            nFilters =nc-2 % on average, 2 sub-units!
            % calculate STA; 
            idx=1:size(X,1);
            spk_tm = idx(double(spks)==1);
            sta=zeros(nc,6); Y=X';
            for ispktm=spk_tm
                if(ispktm>6)
                sta=sta+Y(:,ispktm-5:ispktm);
                end
            end
            sta=sta/length(spk_tm)
            h=figure;
            plot(sta');

            title(sprintf('CellID: %d',cellID));
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/STA_%d.eps',cellID));
end
         
            % maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
            % maskedMov = maskedMov(:,2:end); spks = spks(1:end-1);
            maskedMov = X(1:end-1,:)'; 
            sd = sqrt(diag(maskedMov*maskedMov'/size(maskedMov,2)));
            maskedMov = maskedMov.*repmat(1./sd,[1,size(maskedMov,2)]);
            
            figure;hist(maskedMov(:));title('Filtered Movie')
            maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

            binnedResponsesbigd = double(spks(1,2:end)');
            interval=10;
         
 
 mov_use=maskedMov;
 filteredStimDim=size(mov_use,1);
 nSU = nFilters;
 
 k_est_log=cell(50,1);
 

su_log=zeros(nc,nc);

 for ifit=1:50
 %[fitGMLM,output] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
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

end
%su_log=su_log/50;

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
                if(su_log(icone,jcone)>0)
            plot(datarun.cones.centers(pair,1),-datarun.cones.centers(pair,2),'LineWidth',su_log(icone,jcone));
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
       
s=hgexport('readstyle','coneGMLM');
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/quadGMLM_%d.eps',cellID));
close all

       %% 