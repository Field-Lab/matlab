%% summary plots SU

%% startup
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
%%
ifgmlm=1; % 0 if nnmf

folder = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/';
icellType=2;

%%
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;%

%% 
ifgmlm=1; % 0 if nnmf

folder = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/';
icellType=2;

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;%



datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun =load_sta(datarun)

%%
ifgmlm=1; % 0 if nnmf
folder = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/detailed/';
icellType=10;

WN_datafile = '2015-03-09-2/data031/data031';
WN_datafile_short='2015-03-09-2/data031/data031';
movie_xml = 'BW-8-6-0.48-11111-40x40';
wrong_xml = 'BW-8-6-0.48-22222-40x40';
stim_length=1800;% in seconds
%% Load cellID


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% make figures

for cellID = [datarun.cell_types{icellType}.cell_ids];
    
    close all
    h=figure('Color','w','PaperSize',[42,7],'PaperPosition',[0 0 42 7]);
    iSU=0;
    for nSU=[2,3,4,5,6,7]   
        iSU=iSU+1;
        subplot(1,6,iSU);
       
     
        if(ifgmlm==1)
            load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        %load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
        su_log=su_log;
        %
        mask = totalMaskAccept;
        sta_dim1 = size(mask,1);
        sta_dim2 = size(mask,2);
        indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
        masked_frame = indexedframe(logical(mask));
        
        % make plot
        ilist= 1:length(masked_frame);
        for ix=x_coord
            for iy=y_coord
                i1 =ilist(indexedframe(ix,iy)==masked_frame);
                if(~isempty(i1))
                    
                    i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
                    
                    if(~isempty(i2))
                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[iy,iy+1],'LineWidth',su_log(i1,i2)/4);
                            hold on;
                            text(ix,iy+0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
                    if(~isempty(i2))
                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[iy,iy-1],'LineWidth',su_log(i1,i2)/4);
                            hold on;
                           % text(ix,iy-0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
                    if(~isempty(i2))
                        if(su_log(i1,i2)>0)
                            plot([ix-1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4);
                            hold on;
                         text(ix-0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    
                    i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
                    if(~isempty(i2))
                        if(su_log(i1,i2)>0)
                            plot([ix+1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4);
                            hold on;
                           text(ix+0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                end
                
            end
        end
    axis square    
   title(sprintf('nSU: %d',nSU));
   hold on
    end
  %   xlim([16,32]);ylim([16,32]);

    if(ifgmlm==1)
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/summary/gmlm/gmlm_%d.pdf',cellID));
    else
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/summary/nnmf/nnmf_%d.pdf',cellID))     
    end
    
end


%% Do Spectral clustering

for cellID = [datarun.cell_types{icellType}.cell_ids];
    
    close all
    h=figure('Color','w','PaperSize',[42,7],'PaperPosition',[0 0 42 7]);
    iSU=0;
    for nSU=[2,3,4,5,6,7]   
        iSU=iSU+1;
        subplot(1,6,iSU);
       
     
        if(ifgmlm==1)
            load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        %load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
        su_log=su_log;
        %
        [label] = cluster_spect(su_log,nSU);
        
        
        %
        [col_label]=distinguishable_colors(20);
        
        mask = totalMaskAccept;
        sta_dim1 = size(mask,1);
        sta_dim2 = size(mask,2);
        indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
        masked_frame = indexedframe(logical(mask));
        
        % make plot
        ilist= 1:length(masked_frame);
        for ix=x_coord
            for iy=y_coord
                i1 =ilist(indexedframe(ix,iy)==masked_frame);
               
                if(~isempty(i1))
                  
                    i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
                    
                    if(~isempty(i2))
                        
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end

                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[iy,iy+1],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                            %text(ix,iy+0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[iy,iy-1],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                           % text(ix,iy-0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
                    if(~isempty(i2))
                          if(label(i1)==label(i2))
                        col=col_label(label(i1),:);  
                        else
                        col=col_label(end,:);
                        end
                        if(su_log(i1,i2)>0)
                            plot([ix-1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                         %text(ix-0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    
                    i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            plot([ix+1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                           %text(ix+0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                  
                end
                
            end
        end
        
       
        
        for ix=x_coord
            for iy=y_coord
              i1 =ilist(indexedframe(ix,iy)==masked_frame);
                if(~isempty(i1))
                plot(ix,iy,'*','Color',col_label(label(i1),:));
                end
            end
        end
        
    axis square    
   title(sprintf('nSU: %d',nSU));
   hold on
    end
  %   xlim([16,32]);ylim([16,32]);

    if(ifgmlm==1)
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/ON Parasol/summary_sc/gmlm/gmlm_%d.pdf',cellID));
    else
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/ON Parasol/summary_sc/nnmf/nnmf_%d.pdf',cellID))     
    end
    
end


%% Do Spectral clustering

for cellID = [1531]%[106,1008,842,1066,1232,1531,1981,2596,2767,2371,2806] %[datarun.cell_types{icellType}.cell_ids];
    
    close all
    %h=figure('Color','w','PaperSize',[42,7],'PaperPosition',[0 0 42 7]);
   
    iSU=0;
    for nSU=4%[2,3,4,5,6,7]   
        iSU=iSU+1;
        h=figure('Color','w');
       
     
        if(ifgmlm==1)
            load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        %load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
        su_log=su_log;
        %
        [label] = cluster_spect(su_log,nSU);
        
        
        %
        [col_label]=distinguishable_colors(20);
        
        mask = totalMaskAccept;
        sta_dim1 = size(mask,1);
        sta_dim2 = size(mask,2);
        indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
        masked_frame = indexedframe(logical(mask));
        
        % make plot
        ilist= 1:length(masked_frame);
        for ix=x_coord
            for iy=y_coord
                i1 =ilist(indexedframe(ix,iy)==masked_frame);
               
                if(~isempty(i1))
                  
                    i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
                    
                    if(~isempty(i2))
                        
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end

                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[41,41]-[iy,iy+1],'LineWidth',su_log(i1,i2)/2,'Color',col);
                            hold on;
                            %text(ix,iy+0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            plot([ix,ix],[41,41]-[iy,iy-1],'LineWidth',su_log(i1,i2)/2,'Color',col);
                            hold on;
                           % text(ix,iy-0.5,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
                    if(~isempty(i2))
                          if(label(i1)==label(i2))
                        col=col_label(label(i1),:);  
                        else
                        col=col_label(end,:);
                        end
                        if(su_log(i1,i2)>0)
                            plot([ix-1,ix],[41,41]-[iy,iy],'LineWidth',su_log(i1,i2)/2,'Color',col);
                            hold on;
                        % text(ix-0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    
                    i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            plot([ix+1,ix],[41,41]-[iy,iy],'LineWidth',su_log(i1,i2)/2,'Color',col);
                            hold on;
                           %text(ix+0.5,iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                  
                end
                
            end
        end
        
         
         for ix=x_coord
            for iy=y_coord
                i1 =ilist(indexedframe(ix,iy)==masked_frame);
               
                if(~isempty(i1))
                  
                    i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
                    
                    if(~isempty(i2))
                        
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end

                        if(su_log(i1,i2)>0)
                           % plot([ix,ix],[iy,iy+1],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                            text(ix,41-(iy+0.5),sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            %plot([ix,ix],[iy,iy-1],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                            text(ix,41-(iy-0.5),sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
                    if(~isempty(i2))
                          if(label(i1)==label(i2))
                        col=col_label(label(i1),:);  
                        else
                        col=col_label(end,:);
                        end
                        if(su_log(i1,i2)>0)
                           % plot([ix-1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                         text(ix-0.5,41-iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                    
                    
                    i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
                    if(~isempty(i2))
                            if(label(i1)==label(i2))
                                col=col_label(label(i1),:);
                            else
                                col=col_label(end,:);
                            end
                        if(su_log(i1,i2)>0)
                            %plot([ix+1,ix],[iy,iy],'LineWidth',su_log(i1,i2)/4,'Color',col);
                            hold on;
                           text(ix+0.5,41-iy,sprintf('%d',su_log(i1,i2)));
                            hold on;
                        end
                    end
                  
                end
                
            end
        end
        
        for ix=x_coord
            for iy=y_coord
              i1 =ilist(indexedframe(ix,iy)==masked_frame);
                if(~isempty(i1))
                plot(ix,41-iy,'.','Color',col_label(label(i1),:),'MarkerSize',20);
                end
            end
        end
        set(gca,'xTick',[]);set(gca,'yTick',[]);
    axis equal    
   title(sprintf('nSU: %d',nSU));
   hold on
        xlim([20,27]);ylim([41-28,41-21]);
% 
%     if(ifgmlm==1)
%        
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/detailed_sc_subset/CellID_%d/gmlm/',cellID),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/detailed_sc_subset/CellID_%d/gmlm/',cellID));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/large/detailed_sc_subset/CellID_%d/gmlm/gmlm_%d_su_%d.pdf',cellID,cellID,nSU));
%     else
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_sc_subset/CellID_%d/nnmf/',cellID),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_sc_subset/CellID_%d/nnmf/',cellID));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_sc_subset/CellID_%d/nnmf/nnmf_%d_su_%d.pdf',cellID,cellID,nSU));
%    
%     end

    end
    
end


xxsta =-repelem(sta,1,1);
%     xxsta = xxsta - mean(mean(xxsta((abs(xxsta)<0.2*max(abs(xxsta(:)))))));
B = bwboundaries(abs(xxsta)>0.2*max(abs(xxsta(:))));
hullidx = convhull(B{6}(:,1),B{6}(:,2),'simplify',true);
xxsta = xxsta/max(abs(xxsta(:)));

hold on;
plot(B{6}(hullidx,1),41-B{6}(hullidx,2),'LineWidth',lw);

%% see individual fits


for cellID = 1531%[106,1008,842,1066,1232,1531,1981,2596,2767,2371,2806] %[datarun.cell_types{icellType}.cell_ids];
    
    close all
    %h=figure('Color','w','PaperSize',[42,7],'PaperPosition',[0 0 42 7]);
   
    iSU=0;
    for nSU=4%[2,3,4,5,6,7]   
        iSU=iSU+1;
       
     
        if(ifgmlm==1)
            load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        %load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
        su_log=su_log;
        %
        [label] = cluster_spect(su_log,nSU);
        
        
        %
        [col_label]=distinguishable_colors(20);
        
       
      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);
for ifit=1:50
    close all
    h=figure('Color','w')
    if(ifgmlm==1)
fitGMLM = fitGMLM_log{ifit};
W=zeros(length(masked_frame),nSU);
    else
    W=W_log{ifit};    
    end
for ifilt=1:nSU
subplot(2,2,ifilt);

if(ifgmlm==1)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
else
u_spatial = reshape_vector(W(:,ifilt)/sum(W(:,ifilt)),masked_frame,indexedframe);
end
%subplot(ceil((nSU+1)/2),ceil((nSU+1)/2),ifilt)
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
%title(sprintf('asm Filter: %d',ifilt));
axis square
end

%     if(ifgmlm==1)
%        
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/',cellID,nSU),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/',cellID,nSU));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/gmlm_%d_su_%d_fit_%d.pdf',cellID,nSU,cellID,nSU,ifit));
%     else
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/',cellID,nSU),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/',cellID,nSU));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/nnmf_%d_su_%d_fit_%d.pdf',cellID,nSU,cellID,nSU,ifit));
%    
%     end
end

    end
    
end


%% cellID1531 figure

 cellID = 1531%[106,1008,842,1066,1232,1531,1981,2596,2767,2371,2806] %[datarun.cell_types{icellType}.cell_ids];
    
    close all
    %h=figure('Color','w','PaperSize',[42,7],'PaperPosition',[0 0 42 7]);
    iidx=1:length(datarun.cell_ids);
   sta = -mean(datarun.stas.stas{iidx(datarun.cell_ids==cellID)}(:,:,:,24),3)';
   
    iSU=0;
    nSU=4%[2,3,4,5,6,7]   
        iSU=iSU+1;
       
     
        if(ifgmlm==1)
          load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        %load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
        su_log=su_log;
        %
        [label] = cluster_spect(su_log,nSU);
        
        
        %
        [col_label]=distinguishable_colors(20);
        
       
      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));
x_coord =x_coord(6:end-4);
y_coord = y_coord(5:end-5);

u_spatial_log = zeros(40,40,nSU);
ifit=4%1:50
    lw=1.5;
    close all
    h=figure('Color','w')
    if(ifgmlm==1)
fitGMLM = fitGMLM_log{ifit};
W=zeros(length(masked_frame),nSU);
    else
    W=W_log{ifit};    
    end
    
    subplot(1,nSU+1,1);
    xxsta =-repelem(sta(x_coord,y_coord),20,20);
%     xxsta = xxsta - mean(mean(xxsta((abs(xxsta)<0.2*max(abs(xxsta(:)))))));
    B = bwboundaries(abs(xxsta)>0.2*max(abs(xxsta(:))));
    hullidx = convhull(B{1}(:,1),B{1}(:,2),'simplify',true);
    xxsta = xxsta/max(abs(xxsta(:)));
    imagesc((1-repmat(xxsta,[1,1,3]))/2);
    hold on;
    plot(B{1}(hullidx,2),B{1}(hullidx,1),'LineWidth',lw);
    title('sta');
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);

    colormap gray
      caxis([0,1]);
  % colorbar
    axis image
    
for ifilt=1:nSU
subplot(1,nSU+1,ifilt+1);

if(ifgmlm==1)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
else
u_spatial = reshape_vector(W(:,ifilt)/sum(W(:,ifilt)),masked_frame,indexedframe);
end
%subplot(ceil((nSU+1)/2),ceil((nSU+1)/2),ifilt)
xxu = repelem(u_spatial(x_coord,y_coord),20,20);
xxu = xxu/max(abs(xxu(:)));

imagesc(-xxu);
   hold on;
    plot(B{1}(hullidx,2),B{1}(hullidx,1),'LineWidth',lw);
     caxis([-1,1]);
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
colormap gray
%colorbar
%title(sprintf('asm Filter: %d',ifilt));
axis image
title('Sub-unit')
end
print(h,'-clipboard','-dpdf');
%     if(ifgmlm==1)
%        
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/',cellID,nSU),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/',cellID,nSU));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/SU_%d/gmlm_%d_su_%d_fit_%d.pdf',cellID,nSU,cellID,nSU,ifit));
%     else
% if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/',cellID,nSU),'dir');
%     mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/',cellID,nSU));
% end
%    print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/nnmf/SU_%d/nnmf_%d_su_%d_fit_%d.pdf',cellID,nSU,cellID,nSU,ifit));
%    
%     end



    


