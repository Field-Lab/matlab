path = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/'
load([path,'/Off_type1.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cell_choose_num =8 % [3,7,8];
cellsChoose(cell_choose_num)=1; % [25,26,31,23] or [25]
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

n=40;
xx = repmat([1:n]',[1,n]);
yy = repmat([1:n],[n,1]);
xx_mask = xx(mask);yy_mask = yy(mask);
%%


mov = maskedMovdd(mask,:);
resp = binnedSpikeResponses_coll(cellsChoose,:);

sta = mov*resp';
figure;plotSU2(sta,mask,40,40);

%%
gamma=0.001;
gamma1=1*sqrt((gamma^2)/2);
gamma2 = sqrt(gamma^2 - gamma1^2);
sig = sqrt(var(mov(:)));

Tdiff = zeros(size(mov,1),size(mov,1));
dist_pix = zeros(size(mov,1),size(mov,1));

T =size(resp,2);
for ipix = 1:size(mov,1)
    for jpix = 1:size(mov,1)
    u=zeros(size(mov,1),1);
    u(ipix) = gamma1;
    u(jpix) = -gamma2;
   
    if(ipix~=jpix)
    Tdiff(ipix,jpix) =exp((-norm(u)^2 * sig^2) / 2) * ((exp(u'*mov)*resp' - exp(sig^2*norm(u)^2/2)*sum(resp))/T);
    dist_pix(ipix,jpix) = norm([xx_mask(ipix),yy_mask(ipix)] - [xx_mask(jpix) , yy_mask(jpix)]);
    end
    end
end

% distance - Tdiff scatter
m=size(mov,1);
select_pix = ~logical(eye(m,m));
dl = dist_pix(select_pix);
Tl =abs(Tdiff(select_pix));
figure;plot(Tl,dl,'.');

figure;
subplot(1,3,1);
m=size(mov,1);
bd = diag(ones(m-1,1),1) + diag(ones(m-1,1),-1);
hist(abs(Tdiff(logical(bd))),30);
%hist(Tdiff(~logical(eye(37,37))),30);

Tdiff = abs(Tdiff);

A = double(max(Tdiff(:)) - Tdiff);A = A / max(A(:));A= (A + A')/2;
%A = double(Tdiff < 0.006);
nSU_cluster = 2;
[labels,~]=cluster_spect(A,nSU_cluster); % hierarchical clustering?

subplot(1,3,2);
imagesc(Tdiff(labels,labels));
title('Clustered');

subplot(1,3,3);
plot_SU_clustered(reshape(mask,[40,40]),labels,A*30,nSU_cluster,total_mask_log(:,cellsChoose));


%% variance of T(u) with gamma ? 
ndata=20;
T=size(mov,2);
Tdata=zeros(size(mov,1),size(mov,1),ndata);
Tvar =zeros(size(mov,1),size(mov,1));
gamma_list = [0.1,0.5,1,1.5,2,4,6,8,10];
var_log=[];
igamma = 0;
for gamma=gamma_list;
    igamma = igamma+1;
gamma1=gamma/2;
gamma2 = sqrt(gamma^2 - gamma1^2);
sig = sqrt(var(mov(:)));

for idata = 1:ndata
    movd = mov(:,(idata-1)*T/ndata+1:idata*T/ndata);
    respd=resp(:,(idata-1)*T/ndata+1:idata*T/ndata);
    for ipix =1:size(mov,1)
        for jpix = 1:size(mov,1)
            u=zeros(size(mov,1),1);
            u(ipix) = gamma1;
            u(jpix) = -gamma2;
            
            Tdata(ipix,jpix,idata) = exp(-norm(u)^2 * sig^2/2)*(exp(u'*movd)*respd')/(T/ndata);
        end
    end
     Tvar(ipix,jpix) = var(squeeze(Tdata(ipix,jpix,:)));
   
end
  var_log(igamma) = mean(Tvar(:))
end