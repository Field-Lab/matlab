%% Analyze spatial STA . 



%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/code'));

%% Dataset details
WN_datafile = '2015-03-09-2/d05-27-norefit/data020-from-d05-d27/data020-from-d05-d27';
WN_datafile_short='2015-03-09-2/d05-27-norefit/data020-from-d05-d27/data020-from-d05-d27';
stim_length=30%1799;% in seconds


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_neurons(datarun);

% Load NSEM movie

%% Reference dataset 
WN_datafile = '2015-03-09-2/d05-27-norefit/data018-from-d05-d27/data018-from-d05-d27';
WN_datafile_short='2015-03-09-2/d05-27-norefit/data018-from-d05-d27/data018-from-d05-d27';
stim_length=30%1799;% in seconds


datarun2=load_data(WN_datafile)
datarun2=load_params(datarun2)
datarun2=load_neurons(datarun2);
cellIDs = [datarun2.cell_types{2}.cell_ids,datarun2.cell_types{1}.cell_ids];
%cellIDs = [datarun2.cell_types{1}.cell_ids];
%% Load spikes 

for ref_cell=1:length(cellIDs)
condDuration = stim_length;
cond_str{1}='nsem';
cellID=cellIDs(ref_cell);
neuronPath=datarun.names.rrs_neurons_path;
nConditions=1;
[spkColl,spkCondColl,h]=extract_spikes(cellID,nConditions,condDuration,cond_str,neuronPath)

spks{ref_cell}.spkColl = spkColl;
end

for ref_cell=1:length(cellIDs)
spkColl = spks{ref_cell}.spkColl;
binSz=1/120; len=stim_length/binSz;
spkMat = makeSpikeMat(spkColl,binSz,len);
spks{ref_cell}.spkMat = spkMat;
end
%% Make population activity
imov=1;
binStart = floor(((imov-1)+0.2)/binSz);
binEnd = ceil(((imov-1)+0.8)/binSz);

nTrials = size(spks{ref_cell}.spkMat,1);
clear population
for itrial=1:nTrials
    for icell=1:length(cellIDs)
       population{itrial}.activity(icell,:) = spks{icell}.spkMat(itrial,binStart:binEnd);
    end
end

%% estimate distance
data=[];
for useTrial=1:19
data=[data;population{useTrial}.activity']; % Each row is a data point
end

nCities=size(data,1);
measuredDist=zeros(nCities,nCities);

for icity=1:size(data,1)
    for jcity=1:size(data,1)
   measuredDist(icity,jcity)=norm(data(icity,:)-data(jcity,:));
    end
end


%% Make graph
neighbourSz=20;
E=zeros(nCities);
Edge01 = zeros(nCities);

for icity=1:nCities
distances = measuredDist(icity,:);
[V,idx] = sort(distances,'ascend');
E(icity,idx(1:neighbourSz+1))=V(1:neighbourSz+1);
E(idx(1:neighbourSz+1),icity)=V(1:neighbourSz+1);
Edge01(idx(1:neighbourSz+1),icity)=1;
Edge01(icity,idx(1:neighbourSz+1))=1;
end
Edge01 = logical(Edge01);

%% Calculate minimum distances
% 
% % Dykstra's shortest path algorithm
% idx=1:nCities;
% E_d=zeros(nCities,nCities);
% 
% parfor icity=1:nCities
%     icity
%  E_d(icity,:) = minDistCalc(E,Edge01,icity);  
% end
% 
% E_d(logical(diag(ones(nCities,1))))=0;
% 
% E_d=(E_d+E_d')/2;

%% FastFloyd from internet.
E_d2 = E;
E_d2(Edge01==0)=Inf;
E_d2(logical(diag(ones(size(E,1),1))))=0;
E_d2 = FastFloyd(E_d2);

%E_d2= E_d2(E_d2~=Inf);
E_d2=(E_d2+E_d2')/2; 
% E_d2(E_d2==Inf)=max(E_d2(E_d2~=Inf));
%% 

P=eye(nCities,nCities) - ones(nCities,1)*ones(1,nCities)/nCities;
Q = -P*E_d2*P;
[Eig,lam]=eig(Q);

figure;
plot(diag(abs(lam)),'*')

X = Eig(:,1:2)*abs(lam(1:2,1:2))^(1/2);

col = 'rgbcmyk';

figure;
nTrials=length(population);
trialLen = size(data,1)/nTrials;
for itrial=1:19
% plot(X((itrial-1)*trialLen+1:itrial*trialLen,1),X((itrial-1)*trialLen+1:itrial*trialLen,2),'r*');
% hold on;
%plot(X((itrial-1)*trialLen+1:itrial*trialLen,1),X((itrial-1)*trialLen+1:itrial*trialLen,2),col(itrial));
hold on;
x1=X((itrial-1)*trialLen+1:itrial*trialLen,1)';
x2 = X((itrial-1)*trialLen+1:itrial*trialLen,2)';
z= zeros(trialLen,1)';
col = 1:trialLen; 

surface([x1;x1],[x2;x2],[z;z],[col;col] ,'facecol','no',...
        'edgecol','interp',...
        'linew',2);
end


figure;
nTrials=length(population);
trialLen = size(data,1)/nTrials;
for itrial=1:19
plot(X((itrial-1)*trialLen+1:itrial*trialLen,1),X((itrial-1)*trialLen+1:itrial*trialLen,2),'r*');
hold on;
plot(X((itrial-1)*trialLen+1:itrial*trialLen,1),X((itrial-1)*trialLen+1:itrial*trialLen,2));
hold on;
end

itrial=1;
figure;
for itime=2:trialLen
hold on;
x1=X((itrial-1)*trialLen+1:(itrial-1)*trialLen+itime,1)';
x2 = X((itrial-1)*trialLen+1:(itrial-1)*trialLen+itime,2)';
z= zeros(length(x1),1)';
col = 1:itime; 

surface([x1;x1],[x2;x2],[z;z],[col;col] ,'facecol','no',...
        'edgecol','interp',...
           'linew',2);
       xlim([-5,5]);ylim([-5,5]);
       pause(0.1);
end




%% Dimensionality reduction

data2=data-repmat(mean(data,1),[size(data,1),1]);
C=data2'*data2;
[Eig2,lam2]=eig(C);
[lam_sort,ii]=sort(diag(lam2),'descend');
Eig2=Eig2(:,ii);

figure;
plot(lam_sort,'*');
basis = Eig2(:,1:2)*diag(lam_sort(1:2))^(1/2);

X = data2*basis;