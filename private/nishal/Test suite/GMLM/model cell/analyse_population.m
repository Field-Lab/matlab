load('/Volumes/Lab/Users/bhaishahster/GMLM_fits/ts_population8.mat');

nfits=50;
pairs = zeros(size(K,1),size(K,1));
Bavg = zeros(size(K,1),size(B,2));
for ifit=1:nfits
K = gather(fitASM_pop_log{ifit}.K);
B = gather(fitASM_pop_log{ifit}.B);
su = double(K==repmat(max(K,[],2),[1,size(K,2)]));
pairs = pairs + su*su';

Bavg = Bavg + su*exp(B);
end
pairs =pairs /nfits;
Bavg = Bavg/nfits;


figure('Color','w');
imagesc(pairs);
axis image
colormap gray

figure('Color','w');
imagesc(Bavg');
axis image
colormap gray

lbl = cell(14,1);
for ipix=1:14
lbl{ipix,1} = sprintf('pix %d',ipix);
end

leg=cell(1,Ns);
for isu=1:Ns
   leg{1,isu} = sprintf('extracted su %d',isu);
end

for ifit=1:50

spider((fitASM_pop_log{ifit}.K),'fit',[],lbl,leg);
pause
end