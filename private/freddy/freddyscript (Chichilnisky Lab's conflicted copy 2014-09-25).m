%% Script to look for bidirectional axon activation in 512 stim data - Freddy
%Kellison-Linn

%%
figure; plot(1:1:10, 10:-1:1,'Color',[1 0 0]);
ylim([0 20]); 
set(gca,'YLim',[0 10]); 
title('my title'); 
xlabel('X Label')
% ylabel('Y Label'); 

% Useful plots
figure; imagesc(eye(10)); colorbar; 
axis image; % Makes pixels square
axis xy;    %Changes origin
axis ij; 



axonBundleVolts = zeros(512, 1);
axonBundleVolts = axonBundleVolts';

threshByPattern = NaN(512, 1)




[~, electrodeArray] = ei2matrix(axonBundleVolts)
figure; imagesc(electrodeArray)

cmap = colormap(flipud(gray));
%cmap(1,:) = [0,0,0]

electrodeArray = electrodeArray / 16.05148356;

figure; imagesc(electrodeArray); colorbar; colormap(cmap); title('Approx. # of axons'); axis image;

figure; imagesc(stimThresh512Double); colorbar; colormap(cmap); caxis([0 4]); title('double'); axis image
figure; imagesc(stimThresh512Single); colorbar; colormap(cmap); caxis([0 4]); title('single'); axis image



playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2012-09-24-3/data008/',9,'saveMovie',false)