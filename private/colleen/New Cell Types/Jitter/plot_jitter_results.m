number =545;
load(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026/Cell ', num2str(number), '.mat'])
temp = permute(temp, [2,1,3,4]);
plot_sta_(temp)
title({'2016-02-17-6 data026' ;['Cell ', num2str(number)]})
axis off
