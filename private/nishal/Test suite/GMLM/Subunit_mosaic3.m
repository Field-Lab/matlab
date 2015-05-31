
WN_datafile = '2015-03-09-2/data031/data031';
WN_datafile_short='2015-03-09-2/data031/data01';
movie_xml = 'BW-8-6-0.48-11111-40x40';
stim_length=1800;% in seconds

%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

%%
ref_cell=0;
xx_sum=zeros(400,400);
for cellID =datarun.cell_types{2}.cell_ids(1:end);
ref_cell=ref_cell+1;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Off parasol/Cell%d_full.mat',cellID));

%figure;

for nSU=1:3
    
    fitGMLM_full2 = fitGMLM_full2_log{nSU};
      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);
h=figure;
for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM_full2.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
subplot(2,2,ifilt)
imagesc(repelem(u_spatial(x_coord,y_coord),10,10));
%imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end
  if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Off parasol/CellID_%d',cellID)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Off parasol/CellID_%d',cellID));
  end
  hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Off parasol/CellID_%d/nSU%d.eps',cellID,nSU));
       
% [v,I] = max(u_spatial_log,[],3);
% 
% xx=I.*(v>0.2);
% xxBig = makeBoundarySU(xx,nSU);
% %xx=xx(x_coord,y_coord);
% %figure;
% xx_sum=xx_sum+xxBig;
% 
% xxBig2 = makeBoundarySU(double(xx~=0),nSU);
% xx_sum(xxBig2~=0)=3 ;
% 
% %subplot(2,2,nSU);
% %imagesc(xx);
% %title(sprintf('Num SU : %d',nSU));
% try
% for iSU=1:nSU
% [r,c]=find(xx==iSU);
% k=convhull(r,c,'simplify',true);
% su(ref_cell).ch{iSU}=[r(k),c(k)];
% end
% 
% [r,c]=find(xx>0);
% k=convhull(r,c,'simplify',true);
% su(ref_cell).cellBoundary = [r(k),c(k)];
% catch
% end
end
end

%% 
figure;
for ref_cell=1:31
for iSU=1:nSU
    try
    plot(su(ref_cell).ch{iSU}(:,1),su(ref_cell).ch{iSU}(:,2),'b');
    hold on;
    catch
    end
end
plot(su(ref_cell).cellBoundary(:,1), su(ref_cell).cellBoundary(:,2));
end
