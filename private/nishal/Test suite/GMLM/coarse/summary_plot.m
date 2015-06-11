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

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;%
%% Load cellID


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%% make figures

for cellID = 1531%[datarun.cell_types{icellType}.cell_ids];
    
    close all
    h=figure('Color','w','PaperSize',[45,6],'PaperPosition',[0 0 45 6]);
    iSU=0;
    for nSU=[2,3,4,5,6,7]
        iSU=iSU+1;
        subplot(1,6,iSU);
        if(ifgmlm==1)
            load(strcat(folder,sprintf('CellID_%d/gmlm/Cell%d_gmlm_su_%d.mat',cellID,cellID,nSU)));
        else
            load(strcat(folder,sprintf('CellID_%d/nnmf/Cell%d_nnmf_su_%d.mat',cellID,cellID,nSU)));
        end
        load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
        
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
    end
    xlim([0,40]);ylim([0,40]);
    if(ifgmlm==1)
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/summary/gmlm/gmlm_%d.pdf',cellID));
    else
   print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/summary/nnmf/nnmf_%d.pdf',cellID))     
    end
    
end