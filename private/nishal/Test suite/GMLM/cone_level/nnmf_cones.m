    
      %% Load alex cone data
      addpath(genpath('../../../Test suite/'));
          %  load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');
          % close all
         % load('/Volumes/Analysis/2011-12-13-2/subunits/data008_2/conepreprocess.mat');
            % Good cell - WN , cell ID: 166 thr=0.3
            
           % load('/Volumes/Analysis/2011-12-13-2/subunits/data009/conepreprocess.mat');
           % Good cell - Voronoi, cell ID : 168 ,thr =0.5 , cells=  [3736,1351,2251,3586,4576,5162,1321];
            
          % load('/Volumes/Analysis/2011-12-13-2/subunits/data008/conepreprocess.mat');
            % thr = 0.33
            
            
            load('/Volumes/Analysis/2010-09-24-1/subunits/data036-from-d05_36/conepreprocess.mat');
            thr=0.5; % [827,1157,1907,2341,2851,3346,3961,4081,4306,4460,4637,4681,4862,4876,4936,4998,5011,5176,5191,5251,5476,5701,5851,5911,6076,6151,6241,6332,7036];
            
            
            coneidx=1:size(datarun.cones.weights,1);
            iidx=1:length(datarun.cell_ids);
           for  cellID=[5251,5476,5701,5851,5911,6076,6151,6241,6332,7036];
             close all
             
            ref_cell=iidx(datarun.cell_ids == cellID)%166%1:10;% 4 is interesting.
            %cellID=datarun.cell_ids(ref_cell)
           
            
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
%     binnedResponsesbigd = double(spks(1,2:end)');
%     maskedMov = X(1:end-1,:)';      
    
%     % 3rd from last
%        binnedResponsesbigd = double(spks(1,3:end)');
%     maskedMov = X(1:end-3,:)';   
  %% for voronoi, filter the movie with STA
            binnedResponsesbigd = double(spks(1,1:end)');
            tf=mean(sta(:,end:-1:1),1);
            maskedMov= 40*filterMov_cone(X',logical(ones(nc,1)),squeeze(tf));
            %%
            sd = sqrt(diag(maskedMov*maskedMov'/size(maskedMov,2)));
            maskedMov = maskedMov.*repmat(1./sd,[1,size(maskedMov,2)]);
            
            figure;hist(maskedMov(:));title('Filtered Movie')
            maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

            % Partition the data
            
            trainFrac = 1; % default by Jeremy
            subDivide = 10; % default by Jeremy


%%%% Partition of the data
 [train test] = partitionDat(maskedMov',trainFrac,subDivide,1);
train=train';test=test';

           %%%% Partition of the data
 [resptrain resptest] = partitionDat(binnedResponsesbigd,trainFrac,subDivide,1);

 
            
            %%

 interval=10;
 mov_use=train; 
 resp_use = resptrain;
 
 filteredStimDim=size(mov_use,1);
 for   nFilters = [nc-6: nc-0] % on average, 2 sub-units!
     if(nFilters<1)
     continue;
     end
 nSU = nFilters;

 STX = mov_use(:,resp_use~=0);
 
 STX = STX-min(STX(:));
 STX = STX/max(STX(:));
 
 su_log=zeros(nc,nc);
 W_log=cell(50,1);
 %STX=gpuArray(STX);
 for ifit=1:50
     ifit
     tic;
 [W,H] = nnmf(STX,nSU); % fast as such, but of smaller matrices (~7x20000), gpu is much slower! 
 toc;
 %W=gather(W);
 
 %tol=10^-2;
 %[W,H]=nnmf_gpu(STX,nSU,tol);
 
 k_est=W';
 su=double((k_est)==repmat(max((k_est),[],1),[nSU,1]));
 su_log = su_log + su'*su;
 W_log{ifit}=W;
 end
 
 %%
    icnt=0;
    for isubunits=1:size(su,1)
        if(sum(su(isubunits,:))~=0)
            icnt=icnt+1;
            cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
        end
    end


       
       %
       col='rgbcmykmrgb'
%        figure;
%        for isu=1:length(cells(ref_cell).sub_units)
%            for icone=cells(ref_cell).sub_units(isu).indices
%        plot(datarun.cones.centers(icone,1),-datarun.cones.centers(icone,2),strcat(col(isu),'*'));
%        axis equal
%        hold on
%            end
%            plot(datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,1) ,-datarun.cones.centers(cells(ref_cell).sub_units(isu).indices,2),col(isu));
%            hold on;
%        end

       
       
      h=figure;
      clear icone
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
       
       save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/nnmf/nnmf_%d_su_%d.mat',cellID,nFilters),'W_log','su_log');
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2010-09-24-1/data036/nnmf/nnmf_%d_su_%d.eps',cellID,nFilters));

 end
           end            
%% Cone interaction 
paircones=[4,5]
% plot_cone_interaction(resptrain , mov_use, paircones);
%plot_cone_interaction(binnedResponsesbigd,maskedMov,paircones);
save(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex//2010-09-24-1/data036/cell_%d_interaction_alex.mat',cellID),'binnedResponsesbigd','maskedMov','paircones');
