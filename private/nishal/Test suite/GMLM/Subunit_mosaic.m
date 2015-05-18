
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('~/Nishal/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/private/nishal'));
addpath(genpath('~/Nishal/matlab/code'));
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
%% Dataset details
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% in seconds

%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);
  c = repmat(1:40,[40,1]);
  r = repmat([1:40]',[1,40]);
  Adj_idx = reshape(1:40*40,[40,40]);
  
  cellIDs = datarun.cell_types{2}.cell_ids;

nCells=21;
for ref_cell=1:nCells%1:nCells
    cellID = cellIDs(ref_cell);
%% Significant stixels
        idx=1:length(datarun.cell_ids);
        matlab_cell_ids = idx(datarun.cell_ids==cellID);
        stas=datarun.stas.stas(matlab_cell_ids);
      
        
        % Load STAs
        stas_new=cell(length(stas),1);
        for icell=1:length(stas)
            st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
            for itime=1:size(stas{1},4)
                st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
            end
            stas_new{icell}=st_temp;
        end
        stas=stas_new;
        
        % Used in movie post process
        cell_params.STAlen=14;
        [stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs(stas,cell_params);
        pixr = r(totalMaskAccept2>0);
        pixc = c(totalMaskAccept2>0);
        Adj_idx_list= [];
        for ipix=1:length(pixr)
        Adj_idx_list = [Adj_idx_list;Adj_idx(pixr(ipix),pixc(ipix))];
        end
     %%   
cell_data = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/Cell%d.mat',cellIDs(ref_cell)));
A_log=zeros(3600,3600);
for itrial=1:50
    itrial
   filters = cell_data.fitGMLM_log(itrial).fitGMLM.Linear.filter; 
   mask = cell_data.totalMaskAccept;
   x_coord = cell_data.x_coord;
   y_coord = cell_data.y_coord;
   totalMaskAccept = cell_data.totalMaskAccept;
   sta_dim1 =40;
   sta_dim2 = 40;
   indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
   masked_frame = indexedframe(logical(mask));

   u_spatial_log = zeros(40,40,nSU);

% figure;
for ifilt=1:nSU
subplot(2,2,ifilt)
u_spatial = reshape_vector(filters{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
% imagesc(u_spatial(x_coord,y_coord));
% colormap gray
% colorbar
% title(sprintf('GMLM Filter: %d',ifilt));
% axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

for ipix=1:length(pixr)
    for jpix=1:length(pixr)
        if(I(pixr(ipix),pixc(ipix)) == I(pixr(jpix),pixc(jpix)))
   A_log(Adj_idx(pixr(ipix),pixc(ipix)) ,  Adj_idx(pixr(jpix),pixc(jpix)))=  A_log(Adj_idx(pixr(ipix),pixc(ipix)) ,  Adj_idx(pixr(jpix),pixc(jpix)))+1;
        end
    end
end

end

A_cell_log(ref_cell).A_log=A_log;


 
      h=figure;
      nc=40;
       
       for iconex=1:nc
           for iconey=1:nc
                for jconex=1:nc
                    for jconey=1:nc
                    iidx = Adj_idx(iconex,iconey);
                    jidx = Adj_idx(jconex,jconey);
                    if(A_log(iidx,jidx)>0 & A_log(iidx,jidx)>40)
                    plot([iconex;jconex],[iconey;jconey],'LineWidth',(A_log(iidx,jidx)-40)/4);
                    hold on;
                    end
                    end
                end
           end
       end
       
       
       for iconex=x_coord
           for iconey=y_coord
        scatter(iconex ,iconey,20,'filled','r');
         hold on;
           end
       end
      title(sprintf('Cell ID: %d',cellID));
      hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/SU_grid_Cell%d.eps',cellID)); 
      
      %% SU and Hull
      A_log_slice =A_log(Adj_idx_list,Adj_idx_list);
      
      % Spectral Clustering
      A = A_log_slice;
      npts = size(A,2);
d=sum(A,2);
D=diag(d);
L = eye(npts,npts) - ((inv(D.^(1/2)))*A*(inv(D.^(1/2))));
% L =D-A;
L=(L+L')/2;
L=sparse(L);

[E,lam]=eigs(L,5,'sa');
lam=diag(lam);
[lam_sort,idx]=sort(lam,'ascend');

figure;
plot(lam_sort,'*');

E=E(:,idx);

kf=3%input('Enter number of classes');
U=E(:,1:kf);
for ipt=1:npts
U(ipt,:)=U(ipt,:)/norm(U(ipt,:));
end


[label,C]=kmeans(U,kf);
A_class=zeros(40,40);
for ipix=1:length(pixr)
A_class(pixr(ipix),pixc(ipix))=label(ipix);
end
A_cell_log(ref_cell).A_class = A_class;

A_class_big=zeros(400,400);
for ipix=1:40
    for jpix=1:40
    A_class_big((ipix-1)*10+1 : (ipix-1)*10+10 , (jpix-1)*10+1 : (jpix-1)*10+10)= A_class(ipix,jpix);
    end
end
h=figure;
imagesc(A_class_big);
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/SU_Cell%d.eps',cellID));
close all
end

%% 
A_sum = zeros(3600,3600);
for icell=1:nCells
A_sum = A_sum + A_cell_log(icell).A_log;
end
Adj_idx_list = find(diag(A_sum)>0);

  A_log_slice =A_sum(Adj_idx_list,Adj_idx_list);
      pixr=[];pixc=[];
      for ipix=1:length(Adj_idx_list)
      [r,c]=find(Adj_idx == Adj_idx_list(ipix));
          pixr=[pixr;r];
          pixc=[pixc;c];
      end

      % Spectral Clustering
      A = A_log_slice;
      npts = size(A,2);
d=sum(A,2);
D=diag(d);
L = eye(npts,npts) - ((inv(D).^(1/2))*A*(inv(D).^(1/2)));
% L =D-A;
L=(L+L')/2;
L=sparse(L);

[E,lam]=eigs(L,100,'sa');
lam=diag(lam);
[lam_sort,idx]=sort(lam,'ascend');

figure;
plot(lam_sort,'*');

E=E(:,idx);

kf=input('Enter number of classes');
U=E(:,1:kf);
for ipt=1:npts
U(ipt,:)=U(ipt,:)/norm(U(ipt,:));
end


[label,C]=kmeans(U,kf);
A_class=-1+zeros(40,40);
for ipix=1:length(pixr)
A_class(pixr(ipix),pixc(ipix))=label(ipix);
end
A_cell_log(ref_cell).A_class = A_class;

A_class_big=zeros(400,400);
for ipix=1:40
    for jpix=1:40
    A_class_big((ipix-1)*10+1 : (ipix-1)*10+10 , (jpix-1)*10+1 : (jpix-1)*10+10)= A_class(ipix,jpix);
    end
end
h=figure;
imagesc(A_class_big);
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/Overall_SU.eps',cellID));

%% 

      h=figure;
      nc=40;
       col = 'rgbcmykrgbcmyk';     
       for iconex=1:nc
           for iconey=1:nc
        scatter(iconex ,iconey,10,'filled','r');
         hold on;
           end
       end
       
      for ref_cell=1:13
       A_log = A_cell_log(ref_cell).A_log;
       for iconex=1:nc
           for iconey=1:nc
                for jconex=1:nc
                    for jconey=1:nc
                    iidx = Adj_idx(iconex,iconey);
                    jidx = Adj_idx(jconex,jconey);
                    if(A_log(iidx,jidx)>0 & A_log(iidx,jidx)>40)
                    plot([iconex;jconex],[iconey;jconey],'LineWidth',(A_log(iidx,jidx)-40)/5,'Color',col(ref_cell));
                    hold on;
                    end
                    end
                end
           end
       end
       
 
      end
