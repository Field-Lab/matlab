
%%
path = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/'
load([path,'/Off_type1.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cellsChoose([3,7,8])=1; % [25,26,31,23] or [25]
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;
%%
% mask avaialble

mask = sum(total_mask_log(:,cellsChoose),2)~=0;
Y_C = binnedSpikeResponses_coll(cellsChoose,:);
[r,c] =find(reshape(mask,[40,40])');

figure;
icnt = 0;
for icell=1:3
    for jcell=1:3
        icnt=icnt+1;
subplot(3,3,icnt);
sta{icell,jcell} = reshape(maskedMovdd*(double(Y_C(icell,:)'>0 & Y_C(jcell,:)'>0)),[40,40]);
xx = sta{icell,jcell}';
imagesc(xx(min(r):max(r),min(c):max(c)));
 colormap gray
colorbar
title(sprintf('Cells %d - %d',icell,jcell));
axis image
    end
end
suptitle({'Correlated STA'});

%%

ifit=80;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_Type1_fit_%d.mat',ifit));
[B_use,h]= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));


K = fitASM_pop.K;
B = fitASM_pop.B;
X = maskedMovdd(mask,:);
T = size(X,2);

alpha= cell(3,1);
alpha_max=cell(3,1);
for icell=1:3
cell_inp = exp(K'*X +repmat(B(:,icell),[1,T]));
alpha{icell} = cell_inp ./repmat(sum(cell_inp,1),[Ns,1]);
alpha_max{icell} = double(alpha{icell} == repmat(max(alpha{icell},[],1),[Ns,1]));
end

[[1:Ns]',exp(fitASM_pop.B)]

K=gather(fitASM_pop.K);
K_su = K == repmat(max(K,[],2),[1,size(K,2)]);


% Response curve analysis
cell_a =2;
isu = 3;
jsu=4;
% combined_a = alpha{cell_a}(isu,:) + alpha{cell_a}(jsu,:);
% cond2 =  (combined_a> prctile(combined_a,80));
% commonSU_max = cond2;
commonSU_max = (alpha_max{cell_a}(isu,:)==1) | (alpha_max{cell_a}(jsu,:)==1);
commonSU_max_sta = maskedMovdd*double(commonSU_max)';
[r,c] =find(repelem(reshape(mask,[40,40])',20,20));
xx = (reshape(commonSU_max_sta,[40,40])');

figure
imagesc(repelem(xx,20,20));
colormap gray
hold on;
for icell=1:3
plot(B_use{icell}(:,1),B_use{icell}(:,2),'LineWidth',2);
hold on;
end
xlim([min(c),max(c)]);
ylim([min(r),max(r)]);

% aa = gather(alpha_max{cell_a}(isu,:)');
% (aa(1:end-1)==aa(2:end)) & aa(1:end-1)==1

X_m = X(:,commonSU_max );
figure;
plot(10*mean(X_m,2),'--*');
hold on;
plot(fitASM_pop.K(:,isu),'--*');
hold on;
plot(fitASM_pop.K(:,jsu),'--*');

X_m_sta = mean(X_m,2);

cidx = 1:size(K,1);
figure;
icnt=0;
thrr=0.5;
for ref_pix=cidx(logical(K_su(:,jsu)) & (X_m_sta > thrr*max(X_m_sta)));
    icnt=icnt+1
    subplot(3,4,icnt);
pix1=cidx(logical(K_su(:,isu)) & (X_m_sta > thrr*max(X_m_sta)));
pix2 = cidx(logical(K_su(:,jsu))& (X_m_sta > thrr*max(X_m_sta)));%[75,76,85,86,87];
plot_cone_interaction2(Y_C(2,commonSU_max )',X_m,ref_pix,pix1,pix2,false);
title(sprintf('%d : %d',ref_pix,mean(X_m(ref_pix,:))));
end


figure;
icnt=0;
su = [1,1,2,2];
pix_l =[1,2,1,2]; 
pix1=[46,32];
pix2 = [11,20];%[75,76,85,86,87];

for ref_pix=[pix1,pix2]
icnt=icnt+1;
subplot(2,2,icnt);
plot_cone_interaction2(Y_C(2,commonSU_max )',X_m,ref_pix,pix1,pix2,false);
%title(sprintf('%d : %d',ref_pix,mean(X_m(ref_pix,:))));
title(sprintf('SU %d : Pix %d',su(icnt),pix_l(icnt)));
end
%% overlay with STA from single cone run .. 

%% Multiple fit analysis
K_log = cell(50,1);
B_log = cell(50,1);
K_su_log = cell(50,1);
d=99;
K_su_add = zeros(d,d);

icnt=0;
for ifit = 22:54
    icnt=icnt+1;
data = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_type1_fit_%d.mat',ifit));
K_log{icnt} = gather(data.fitASM_pop.K);
B_log{icnt} = gather(data.fitASM_pop.B);
K_su = double(K_log{icnt} == repmat(max(K_log{icnt},[],2),[1,size(K_log{icnt},2)]));
K_su_log{icnt} = K_su;
K_su_add = K_su_add + K_su*K_su';
end
K_su_add=K_su_add;
nSU_cluster=15%data.Ns;
[labels,h]=cluster_spect(K_su_add,nSU_cluster); % hierarchical clustering?

figure;
imagesc(K_su_add(labels,labels));
title('Clustered');

figure;
plot_SU_clustered(reshape(data.mask,[40,40]),labels,K_su_add*50/33,nSU_cluster,total_mask_log(:,cellsChoose));


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