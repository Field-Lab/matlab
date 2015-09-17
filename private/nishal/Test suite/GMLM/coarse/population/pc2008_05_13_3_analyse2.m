% pc 2008_05_13_3 individual fits analysis


ifit=87;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_Type1_fit_%d.mat',ifit));

K_log = cell(50,1);
K_su_log = cell(50,1);
K_su_add = zeros(sum(mask_Cell),sum(mask_Cell));

for ifit =1 :50
filters = fitASM_log{ifit}.Linear.filter;
K=zeros(sum(mask_Cell),1);
for isu = 1:Ns
    K(:,isu) = filters{isu};
end

K_log{ifit}=K;
K_su = double(K_log{ifit} == repmat(max(K_log{ifit},[],2),[1,size(K_log{ifit},2)]));
K_su_log{ifit} = K_su;
K_su_add = K_su_add + K_su*K_su';
end


K_su_add=K_su_add;
nSU_cluster=4%data.Ns;
[labels,h]=cluster_spect(K_su_add,nSU_cluster); % hierarchical clustering?

figure;
imagesc(K_su_add(labels,labels));
title('Clustered');

figure;
plot_SU_clustered(reshape(mask_Cell,[40,40]),labels,K_su_add,nSU_cluster,total_mask_log(:,cellsChoose));

%%
% plot max value in su v/s min value in su for each pixel.
K_mm = zeros(size(K_su,1),icnt,2);
Ns = size(K_log{1},2);
K_orth = zeros(Ns,Ns);
d = size(K_log{1},1);
Klog_mat = zeros(33,d,Ns);
for icnt=1:33
for ipix=1:size(K_mm,1)
    K = K_log{icnt}   ;
    
    % max and second largest pixel value
    K_sort = sort(K(ipix,:),'descend');
    K_mm(ipix,icnt,1) = K_sort(1);
    K_mm(ipix,icnt,2) = K_sort(2);
    
   
    % orthogonality of SU
    K_norm = K./repmat(sqrt(sum(K.^2,1)),[d,1]);
    K_orth = K_orth + K_norm'*K_norm;
    
    % positivity
    Klog_mat(icnt,:,:)  = K_norm;% K_log{icnt};
end
end
K_orth = K_orth/33;

figure;
x = K_mm(:,:,1);
y = K_mm(:,:,2);
plot(x(x>0.3*(max(x(:)))),y(x>0.3*max(x(:))),'.')
hold on;
plot([0,2],[0,2]);

% inner product between sub-units, 1 fit
figure;
imagesc(repelem(real(acosd(K_norm'*K_norm)),20,20));colormap gray
set(gca,'xTick',[]);
set(gca,'yTick',[]);
axis image

figure; 
hist(Klog_mat(:),40);set(gca,'yTick',[]);
% plot(N,log(X));