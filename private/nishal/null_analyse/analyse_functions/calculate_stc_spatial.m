function [STA_sp,STC_spatial,u_coll,data] = calculate_stc_spatial(spksGen,mov,STA,mask,time_c)
%% Get Mask
global maskedMov binnedResponses masksz

STAs = zeros(size(STA,1),size(STA,2),3,size(STA,3));
STAs(:,:,1,:) = STA(:,:,end:-1:1);
STAs(:,:,2,:) = STA(:,:,end:-1:1);
STAs(:,:,3,:) = STA(:,:,end:-1:1);

stas{1}=STAs;
cell_params.STAlen= size(STA,3);

if isempty(mask)
[new_stas,totalMaskAccept,CellMasks]=clipSTAs(stas,cell_params);
mask = CellMasks{1};
else
    new_stas{1} = STAs.*repmat(mask,[1,1,3,size(STA,3)]);
CellMasks{1}=mask;    
end
%% Debug
% mask = 0*mask ; 
% mask(25,13)=1;
% mask(26,13)=1;
% mask(27,13)=1;
% mask(25,14)=1;
% mask(26,14)=1;
% mask(27,14)=1;
% mask(24,13)=1;
% 
% mask(25,12)=1;
% mask(26,12)=1;
% mask(27,12)=1;
% 
% mask(25,11)=1;
% mask(26,11)=1;
% mask(27,11)=1;
% 
% mask=mask';
% CellMasks{1}=mask;
%% 
if(isempty(time_c))
time_c = squeeze(sum(sum(sum(new_stas{1},1),2),3));

time_c=(time_c(end:-1:1))/sum(mask(:));
end
%time_c(15:end)=0;

figure;
plot(time_c);

mov=squeeze(mean(mov,3));
Filtlen=size(STA,3);
filtMov=zeros(size(mov,1),size(mov,2),size(mov,3)-Filtlen+1);
    for x=1:size(mov,1)
        for y=1:size(mov,2)
        filtMov(x,y,:)=conv(squeeze(mov(x,y,:)),time_c,'valid');
        end
    end


%% Filter movie using mask 

maskedMov=zeros(sum(CellMasks{1}(:)),size(filtMov,3));

    for itime=1:size(filtMov,3)
     xx=filtMov(:,:,itime);
     maskedMov(:,itime) = xx(logical(CellMasks{1}));
    end

    
maskedSTA=zeros(sum(CellMasks{1}(:)),size(STA,3));

    for itime=1:size(STA,3)
     xx=STA(:,:,itime);
     maskedSTA(:,itime) = xx(logical(CellMasks{1}));
    end
%% Calculate spatial STA . Just sanity check and safe calculations.


Filtdim1=size(maskedMov,1);
binnedResponses = spksGen;
Len = length(spksGen);

% My own STA code 

% STA_sp=zeros(Filtdim1,1);
% idx = 1:Len;
% spkbin=idx(binnedResponses==1);
% 
% for ibin=spkbin
%  
%    if(ibin<size(maskedMov,2))
%   if(mod(ibin,1000)==1)
%        ibin
%   end
% STA_sp=STA_sp+maskedMov(:,ibin)*binnedResponses(ibin);
%    end
% end
% 
% STA_sp=STA_sp/sum(binnedResponses);


%STA_sp from one of the ways I did in spatial nulling stuff?  
STA_sp=zeros(size(maskedSTA,1),1);
for idim=1:size(maskedSTA,1)
STA_sp(idim) = maskedSTA(idim,:)*time_c(end:-1:1);
end

figure;
plot(STA_sp);
      
%% 
xxreSTA =STA_sp;
xxreSTA=xxreSTA/norm(xxreSTA(:));

dim1 = size(xxreSTA,1);
binnedResponsesTrial=binnedResponses;
framesValid = find(binnedResponsesTrial>0);

% flatten movie

reSTC=zeros(dim1,dim1);
for iframe=framesValid'
    iframe
    if(iframe<=size(maskedMov,2))
  xx= maskedMov(:,iframe);

 % xxnew = (xx-((xx'*xxreSTA)*xxreSTA));
  xxnew=xx; 
 display('STC screwed up')
reSTC=reSTC+ xxnew*xxnew'*binnedResponsesTrial(iframe);
    end
end


reSTC = reSTC / (sum(binnedResponses)-1);

STC_spatial=reSTC;

%%  SVD for STC_spatial

[V,D]=eig(STC_spatial);
[s,I]=sort(abs(diag(D)),'descend');
u=V(:,I);

figure;
plot(s,'*')

% Reconstruct u into RFs
sta_dim1=size(STA,1);
sta_dim2=size(STA,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(CellMasks{1}));

figure; 
subplot(3,2,1);
u_spatial = reshape_vector(STA_sp,masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('STA'));


for ifilt=1:3
subplot(3,2,ifilt+1)
u_spatial = reshape_vector(u(:,ifilt),masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('Filter: %d',ifilt));
end

u_coll=[];
for ifilt=1:size(u,2)
u_spatial = reshape_vector(u(:,ifilt),masked_frame,indexedframe);
u_coll=[u_coll,u_spatial(:)];
end

%% STC -> Similarity matrix.

pos_STC=abs(STC_spatial);
[V,D]=eig(pos_STC);
[s,I]=sort(abs(diag(D)),'descend');
u=V(:,I);

figure;
plot(s,'*')
title('Positive STC');

% Reconstruct u into RFs
sta_dim1=size(STA,1);
sta_dim2=size(STA,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(CellMasks{1}));

figure; 
subplot(3,2,1);
u_spatial = reshape_vector(STA_sp,masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('STA'));


for ifilt=1:4
subplot(3,2,ifilt+1)
u_spatial = reshape_vector(u(:,ifilt),masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('Postitive STC Filter: %d',ifilt));

end

%% Spectral clustering
A = pos_STC;
d = sum(A,2);
D=diag(d);
dim = size(pos_STC,1);
L = eye(dim,dim) - (D^(-0.5))*A*(D^(-0.5));

[V,D]=eig(L);
[s,I]=sort(abs(diag(D)),'ascend');
u=V(:,I);
k=4;
uk = u(:,1:k);

figure;
plot(s,'*');
title('Spectral clustering - laplacian');

for in=1:size(u,1)
uk(in,:) = uk(in,:)/norm(uk(in,:));
end

[idx,C]=kmeans(uk,k);



figure; 
subplot(3,2,1);
u_spatial = reshape_vector(STA_sp,masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('STA'));


for ifilt=1:4
subplot(3,2,ifilt+1)
u_spatial = reshape_vector(double(idx==ifilt),masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('Spectral clustering Filter: %d',ifilt));

end

%% CVX
% 
% data_short=maskedMov';
% dt=1/120;
% resp=binnedResponses;
% n=length(resp);
% 
% tic;
% cvx_begin
% variable A(7,7) semidefinite
% variables b(7,1) c
% subject to 
% 
% lam = data_short*b + c + sum((A*data_short').*data_short',1)';
% lam=lam';
% f=+(sum(exp(0.15*lam)*dt) - sum(0.15*lam(resp==1))) + sum(gamma*abs(A(:)));
% 
% minimize (f)
% cvx_end
% toc;
% 
%  data.A=A;
%  data.b=b;
%  data.c=c;

 %% Optimization using gradient descent - non convex.
%  
% masksz = sum(mask(:));
%  nargs = 1 + masksz + masksz^2 ;
% % x0=rand(nargs,1)/10;
% x0=ones(nargs,1)/10; % initialization very important !! 
% %x0(1) = 0.4481;
% %x0(2:2+masksz-1)=0.8009;
% 
%  options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
% [x,fval] = fminunc(@quad_GLM_B,x0,options);
% 
% c=x(1);
% b=x(2:2+masksz-1);
% B=reshape(x(2+masksz:end),[masksz,masksz]);
% 
% 
%  data.A=B*B';
%  data.b=b;
%  data.c=c;
  
%% Optimization using gradient descent - non convex and fixed b
masksz = sum(mask(:));
global rankB
rankB=4;
 nargs = 1 + 1 + masksz*rankB ;
% x0 = randn(nargs,1);
 %x0=2*((rand(nargs,1)>0.5)-1)/10;
x0=ones(nargs,1)/10+rand(nargs,1)/100; % initialization very important !! 
%x0(1) = 0.4481;
%x0(2:2+masksz-1)=0.8009;
global b

b=[20,20,20,20,20,20,20]';
options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
[x,fval] = fminunc(@quad_GLM_B_b_fixed,x0,options); % could fix one column of B to ones.

c=x(1);
b_sc=x(2);
B=reshape(x(3:end),[masksz,rankB]);


 data.A=B*B';
 data.b=b*b_sc;
 data.c=c;

 %% Optimization using gradient descent - convex.
% masksz = sum(mask(:));
%  nargs = 1 + masksz + masksz^2 ;
% % x0=rand(nargs,1)/10;
% x0=ones(nargs,1)/10;
% %x0(1) = 0.4481;
% %x0(2:2+masksz-1)=0.8009;
% global gamma
% gamma=500;
% options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
% [x,fval] = fminunc(@quad_GLM_A,x0,options);
% 
% c=x(1);
% b=x(2:2+masksz-1);
% A=reshape(x(2+masksz:end),[masksz,masksz]);
% 
% 
%  data.A=A;
%  data.b=b;
%  data.c=c;
%  
end



function u_spatial = reshape_vector(u,masked_frame,indexed_frame)
dim1 = size(indexed_frame,1);
dim2 = size(indexed_frame,2);
u_spatial= zeros(dim1,dim2);

for ielem=1:length(u)
[r,c]=find(masked_frame(ielem)==indexed_frame);
u_spatial(r,c) = u(ielem);
end


end
