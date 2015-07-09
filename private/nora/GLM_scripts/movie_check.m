clear
disp('loading orig movie')
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat')
%block_check_index = randperm(59);

%%
disp('selecting check blocks')
block_check_index = [20 40 52];
for i = 1:3
    blocks{i} = NSEMmovie.fitmovie.movie_byblock{block_check_index(i)}.matrix;
end
clear NSEMmovie

%%
disp('loading new blocks to compare')
for i =1:3
    load(['/Volumes/Lab/Users/Nora/NSEM_Movies/eye-120-3_0-3600/movieblock' num2str(block_check_index(i)) '.mat'])
    blocks_new{i} = movie.matrix;
    clear movie
end

%%
disp('plotting!')
for i =1:3
    for j =7100:7150
        subplot(1,2,1)
        imagesc(blocks{i}(:,:,j));
        axis image
        colormap gray
        subplot(1,2,2)
        imagesc(squeeze(blocks_new{i}(j,:,:)));
        axis image
        colormap gray
        pause(0.01)
    end
end

%%
clear
disp('Loading orig movie')
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat')
disp('loading new movie')
load('/Volumes/Lab/Users/Nora/NSEM_Movies/eye-120-3_0-3600/testmovie.mat')


%%
disp('plotting!')
for j =3500:3600
    subplot(1,2,1)
    imagesc(testmovie.matrix(:,:,j));
    axis image
    colormap gray
    subplot(1,2,2)
    imagesc(squeeze(movie.matrix(j,:,:)));
    axis image
    colormap gray
    pause(0.01)
end

