% Null space simulation suite ! 

%% Generate RFs - subunit model
nSubunits = 4;
Filtdim1 = 6;
Filtdim2 = 6;
Filtlen = 30;

subunits=cell(nSubunits,1);

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,160-158,:,:)=1;
k(323-318,162-158,:,:)=1;
k(323-318,163-158,:,:)=1;
subunit_scale=1%1;
subunits{1}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,159-158,:,:)=1;
subunit_scale=1%1.2;
subunits{2}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,161-158,:,:)=1;
subunit_scale=1%0.5;
subunits{3}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,160-158,:,:)=1;
k(319-318,159-158,:,:)=1;
subunit_scale=1%1.5;
subunits{4}=k*subunit_scale;


figure;
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,1,4));
colormap gray
title(sprintf('Subunit: %d',isubunit));
end

figure;
mask=double((subunits{1}(:,:,1,4)+subunits{2}(:,:,1,4)+subunits{3}(:,:,1,4)+subunits{4}(:,:,1,4))~=0);
imagesc(mask);
% Temporal properties?

scale_one=1;
scale_two=0.25;
tau_one=4;
tau_two=10;
n_filters=6;
t=[0:29];
tf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));
figure;
plot(tf)

tf2=zeros(1,1,1,Filtlen);
tf2(1,1,1,:)=tf;
tf=tf2;
clear tf2
tfRep=repmat(tf,[Filtdim1,Filtdim2,1,1]);
title('Temporal Filter');

subunits{1}=subunits{1}.*tfRep;
subunits{2}=subunits{2}.*tfRep;
subunits{3}=subunits{3}.*tfRep;
subunits{4}=subunits{4}.*tfRep;


% sub-unit weights
subunitWeights=zeros(nSubunits,1);
subunitWeights(1)=1%1;
subunitWeights(2)=1%1.2;
subunitWeights(3)=1%0.7;
subunitWeights(4)=1%1.5;

% sub-unit non-linearity

f= @(x) double(x>0).*(1.6*x);
%f=@(x) exp(0.6*x);
% Ganglion cell non-linearity
N= @(x) exp(0.15*(x));
%N=@(x) double(x>0)*(0.02).*(3.4*x).^2; % use it!
%N=@(x) x-min(x(:));
%N = @(x) 15./(1+exp(-1.5*(x-5)));
% 
% 
% figure;
% for itime=1:30
%     itime
% for isubunit=1:nSubunits
% subplot(2,2,isubunit);
% imagesc(subunits{isubunit}(:,:,1,itime));
% caxis([min(subunits{isubunit}(:)),max(subunits{isubunit}(:))]);
% colormap gray
% colorbar
% %title(sprintf('Subunit: %d',isubunit));
% 
% end
% % pause
% end
%% make WN movie
make_WN_mov

%% generate response
generate_resp


%% Quadratic approximation to sub-unit non-linearity
input = cell_resp(:,2);
X = [ones(length(input),1),input,input.^2];
y=f(input);
a=X\(X'\(X'*y));

figure;
[x,NN]=hist(input,100);
plot(NN,(x/sum(x)) * 500,'b');
hold on;
plot(NN,f(NN),'r');
hold on;
plot(NN,[ones(length(NN),1),NN',NN'.^2]*a,'k');
legend('Input data histogram','sub-unit Nonlinearity','Quadratic approximation');

f_old=f;
f_new=@(x) (a(1)+a(2)*x+a(3)*x.^2);
f=f_new;

%% Generate WN movie
make_WN_mov
%% Generate response 
generate_resp
%% STA/STC calculations
idx=1:length(binnedResponses);
Filtlen=30;
ifSTC=0;
spksGen = binnedResponses;
mov3=zeros(6,6,3,size(mov,3)+Filtlen-1);
mov3(:,:,1,:)=mov2;
mov3(:,:,2,:)=mov2;
mov3(:,:,3,:)=mov2;

mov4=zeros(6,6,3,size(mov,3));
mov4(:,:,1,:)=mov;
mov4(:,:,2,:)=mov;
mov4(:,:,3,:)=mov;

[STA,STC] = calculate_sta_stc(spksGen,mov4,Filtlen,ifSTC);
 
[STA_spatial,STC_spatial,u_coll,data] = calculate_stc_spatial(spksGen,mov3,STA,mask,squeeze(tf));

%% Create null stimulus for the STA

make_WN_mov
weights = [2,2,2,2,2,2,2]%data.b;
A=weights;
[u,s,v]=svd(A,'econ');
Ainv=v*(s^-1)*u';

indexedframe = reshape(1:Filtdim1*Filtdim2,[Filtdim1,Filtdim2]);
masked_frame = indexedframe(logical(mask));
mov_modify=0*mov;
tic;
for iframe=1:size(mov,3)
    if(rem(iframe,1000)==1)
    iframe
    end
    mov_fr=mov(:,:,iframe);
mov_fr=mov_fr(logical(mask));
mov_null=mov_fr-Ainv*A*mov_fr;
mov_null = reshape_vector2(mov_null,masked_frame,indexedframe);
mov_modify(:,:,iframe)=mov_null;
end
toc;

mov=mov_modify;
%% 
% generate response
generate_resp


%% STA/STC calculations
idx=1:length(binnedResponses);
Filtlen=30;
ifSTC=0;
spksGen = binnedResponses;
mov3=zeros(6,6,3,size(mov,3)+Filtlen-1);
mov3(:,:,1,:)=mov2;
mov3(:,:,2,:)=mov2;
mov3(:,:,3,:)=mov2;

mov4=zeros(6,6,3,size(mov,3));
mov4(:,:,1,:)=mov;
mov4(:,:,2,:)=mov;
mov4(:,:,3,:)=mov;

[STA,STC] = calculate_sta_stc(spksGen,mov4,Filtlen,ifSTC);
 
[STA_spatial,STC_spatial,u_coll,data] = calculate_stc_spatial(spksGen,mov3,STA,mask,squeeze(tf));
