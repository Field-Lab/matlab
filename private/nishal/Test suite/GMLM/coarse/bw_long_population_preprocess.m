% Use this code to save things to .mat file to use it in python later.

%load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_Unknown_6902.mat');

% load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_OFFPar_901.mat');
% load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat');

load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/BW-8-1-0.48-11111_RNG_16807/fitmovie_8pix_Identity_8pix.mat');

movie_frames = BWmovie.fitmovie.movie_byblock;

nmovies=length(movie_frames);
movie_full=zeros(80,40,3600*nmovies);
for imov=1:nmovies
    imov
movie_full(:,:,(imov-1)*3600+1:imov*3600) = double(movie_frames{imov}.matrix) - 0.5;
end


%%
cell_list = [467,497,1292,1321,1471,1771,1786,2211,2506,2851,3361,3676,3843,4336,4367,4576,5086,5161,5386,5851,5986,6196,6256,7066,7292,7606]
binnedSpikeResponses_coll = zeros(length(cell_list),size(movie_full,3));
ttf = zeros(30,length(cell_list));
total_mask_log = zeros(80*40,length(cell_list));
for icell_list=1:length(cell_list)
    icell_list
load(sprintf('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/WN_mapPRJ/organizedspikes_OFFPar_%d.mat',cell_list(icell_list)));
spikes_by_block = organizedspikes.block.t_sp_withinblock(2:2:end);

binnedSpikeResponses =[];
for imov=1:nmovies 
    ss{1}=spikes_by_block{imov}*20000; % Remove spikes which are too soon after change of movie.
    ss{1}(ss{1}<20000*30/120)=-10;
    spkMat = makeSpikeMat(ss,1/120,3600);
    binnedSpikeResponses = [binnedSpikeResponses,spkMat];
end

binnedSpikeResponses_coll(icell_list,:)=binnedSpikeResponses;
%% STA

times = 1:length(binnedSpikeResponses);
t_sp=times(binnedSpikeResponses~=0);

STA=zeros(80,40,30);
for it_sp=t_sp(t_sp>40)
STA = STA+movie_full(:,:,it_sp:-1:it_sp-30+1);
end

STA =STA/numel(t_sp);

average_mov = repmat(mean(movie_full,3),[1,1,30]);


STA2 = STA - average_mov;
STA3 = STA./average_mov;

figure;
subplot(4,1,1);
imagesc(STA(:,:,6)');
colormap gray
axis image

subplot(4,1,2);
imagesc(STA2(:,:,6)');
colormap gray;
axis image

subplot(4,1,3);
imagesc(STA3(:,:,6)');
colormap gray;
axis image


subplot(4,1,4);
imagesc(-STA3(:,:,6)');
colormap gray;
axis image

%%
%% Find relevant pixels and tf using NSEM
        stas{1}=zeros(80,40,3,30);
        stas{1}(:,:,1,:) = STA;
        stas{1}(:,:,2,:) = STA;
        stas{1}(:,:,3,:) = STA;
        
         cell_params.STAlen=30;
         cell_params.thres=1.5;
       %[stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs(stas,cell_params);
        [stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs_largestblob(stas,cell_params);
        % fill in the bad entries in totalMaskAccept2
       total_mask_log(:,icell_list) = totalMaskAccept2(:);
%         
       ttf(:,icell_list)=squeeze(mean(mean(STA2.*repmat(totalMaskAccept2,[1,1,30]),1),2));
%       ttf=1000*ttf;
%       %tf=tf.*double(idx<15)';
%      
       figure;
       plot(ttf);
end

ttf = mean(ttf,2);


%% Filter with mask and tf
      mov = movie_full; 
      maskedMovdd= filterMov(mov,logical(0*totalMaskAccept2+1),squeeze(ttf));
      
      save('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/Off_par2.mat','-v7.3','maskedMovdd','binnedSpikeResponses_coll','ttf','total_mask_log');
     