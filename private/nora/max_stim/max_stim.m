load('/Users/Nora/Dropbox/DOVESstimuli_MHT/DOVESfixations_032114.mat')
load('/Users/Nora/Dropbox/DOVESstimuli_MHT/DOVESframes.mat')

%%
imagesc(DOVESframes{1,1}(:,:,1));
axis image
colormap gray

%%
for i = 1:5
    for j = 1:5
        hist(DOVESframes{i,j}(:));
        pause()
    end
end
