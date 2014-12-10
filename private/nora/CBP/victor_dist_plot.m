files=dir('CBPmat/distance_mat/*.mat');

for i=1:length(files)
    load([ 'CBPmat/distance_mat/' files(i).name])
end

centers=0:3:150;

[counts_CBP_rec, ~]=hist(CBP_dist_rec(:),centers);
length_CBP_rec=length(find(CBP_dist_rec(:)));
[counts_CBP_sim, ~]=hist(CBP_dist_sim(:),centers);
length_CBP_sim=length(find(CBP_dist_sim(:)));
[counts_cluster_sim, ~]=hist(cluster_dist_sim(:),centers);
[counts_cluster_rec, ~]=hist(cluster_dist_rec(:),centers);
[counts_cross_sim, ~]=hist(cross_dist_sim(:),centers);
[counts_cross_rec, ~]=hist(cross_dist_rec(:),centers);


figure;
subplot(2,1,1)
plot(centers(2:end), counts_CBP_rec(2:end)/length(find(CBP_dist_rec(:))),'r');
hold on
plot(centers(2:end),counts_cluster_rec(2:end)/length(find(cluster_dist_rec(:))));
plot(centers(2:end),counts_cross_rec(2:end)/length(find(cross_dist_rec(:))),'k');
title('Distribution of differences between recorded spike trains')
legend('CBP','Clustering','Cross')
ylim([0 0.3])

subplot(2,1,2)
hold on
plot(centers(2:end),counts_CBP_sim(2:end)/length(find(CBP_dist_sim(:))),'r');
plot(centers(2:end),counts_cluster_sim(2:end)/length(find(cluster_dist_sim(:))));
plot(centers(2:end),counts_cross_sim(2:end)/length(find(cross_dist_sim(:))),'k');
title('Distribution of differences between simulated spike trains')
legend('CBP','Clustering','Cross')
ylim([0 0.3])

   
