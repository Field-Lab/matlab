       %% Load alex cone data
            load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');


            coneidx=1:size(datarun.cones.weights,1);

            ref_cell=4% 4 is interesting.
            nFilters =5;
            cellID=datarun.cell_ids(ref_cell)
           
            
            spks=datarun.spike_rate(ref_cell,:);
            X =double(datarun.cone_inputs(:,logical(datarun.cones.weights(:,ref_cell))));
            cones = coneidx (logical(datarun.cones.weights(:,ref_cell)));
            nc=length(cones);
            
            % calculate STA; 
            idx=1:size(X,1);
            spk_tm = idx(double(spks)==1);
            sta=zeros(nc,6); Y=X';
            for ispktm=spk_tm
                sta=sta+Y(:,ispktm-5:ispktm);
            end
            sta=sta/length(spk_tm);
         
            % maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
            % maskedMov = maskedMov(:,2:end); spks = spks(1:end-1);
            maskedMov = X(1:end-1,:)'; 
            
            figure;hist(maskedMov(:));title('Filtered Movie')
            maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

            binnedResponsesbigd = double(spks(1,2:end)');
            interval=10;
         
 
 mov_use=maskedMov;
 filteredStimDim=size(mov_use,1);
 nSU = nFilters;
 [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
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

        

su=double((-k_est)==repmat(max(-(k_est)),[nFilters,1]));
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
       end
