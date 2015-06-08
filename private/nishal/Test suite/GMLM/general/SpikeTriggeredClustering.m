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

% %% Dataset details
% WN_datafile = '2014-11-24-3/streamed/data008/data008';
% WN_datafile_short='2014-11-24-3/streamed/data008/data008';
% movie_xml = 'BW-1-6-0.48-11111';
% stim_length=10%1799;% in seconds
% cellID = 301;
% 
% 
% %% Load stimuli
% 
% 
% datarun=load_data(WN_datafile)
% datarun=load_params(datarun)
% datarun=load_sta(datarun);
% datarun=load_neurons(datarun);
% 
% 
% %
% stim_description=movie_xml;
% xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
% dashes=find(stim_description=='-');
% StimulusPars.type=stim_description(1:dashes(1)-1);
% StimulusPars.pixelsize = str2double(stim_description(dashes(1)+1:dashes(2)-1));
% StimulusPars.refreshrate = str2double(stim_description(dashes(2)+1:dashes(3)-1));
% StimulusPars.RNG = str2double(stim_description(dashes(3)+1:dashes(4)-1));
% % try
% %     StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
% %     StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
% % catch
% %     StimulusPars.height = 32; 
% %     StimulusPars.width = 64;
% % end
% StimulusPars.tstim = 1/120;
% fitframes=stim_length*120; % seconds * 120 frames per second / interval
% 
% 
% disp('Loading Stimulus Movies')
% [temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers(datarun.triggers<10), fitframes/StimulusPars.refreshrate);
% temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
% fitmovie_color=zeros(width,height,3,fitframes);
% 
% try
%     StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
%     StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
% catch
%     StimulusPars.height =height; 
%     StimulusPars.width = width;
% end
% 
% 
% for i=1:fitframes
%     fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/StimulusPars.refreshrate));
% end
% 
% clear temp_fitmovie i 
% 
% %% Load spikes 
%     master_idx         = find(datarun.cell_ids == cellID);
%     
%       
%         % Spike loading
%         spikes=datarun.spikes{master_idx};
%     
%         % make STA 3D ? 
%         %glm_cellinfo.WN_STA = squeeze(sum(glm_cellinfo.WN_STA,3)); % Doubt!!!!!!!
%         clear cell_savename
%         
%         % Align the spikes and the movies;
%         spikes_adj=spikes;
%         n_block=0;
%         for i=1:(length(datarun.triggers)-1)
%             actual_t_start=datarun.triggers(i);
%             supposed_t_start=n_block*100/120;
%             idx1=spikes > actual_t_start;
%             idx2=spikes < datarun.triggers(i+1);
%             spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
%             n_block=n_block+1;
%         end
%         clear spikes
%         spike.home=spikes_adj;
%         clear spikes_adj;
% %% 
%  spksGen = zeros(stim_length*120,1);
%         for ispike=1:length(spike.home)
%         spksGen(floor(spike.home(ispike)*120)+1)=1;
%         end
%         mov = fitmovie_color-0.5;
%         Filtlen=30;
        
        %% Load alex cone data
        load('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/conepreprocess.mat');
        coneidx=1:size(datarun.cones.weights,1);
        
       for  ref_cell=1% 4 is interesting.
        
        cellID=datarun.cell_ids(ref_cell);
        spks=datarun.spike_rate(ref_cell,:);
        X =double(datarun.cone_inputs(:,logical(datarun.cones.weights(:,ref_cell))));
        cones = coneidx (logical(datarun.cones.weights(:,ref_cell)));
        
        idx=1:length(spks);
        spk_frames = idx(spks>0);
        STx= X(spk_frames-1,:); %DT 
        STx = STx - repmat(mean(STx,1),[size(STx,1),1]);
        npts = size(STx,1);
A=zeros(npts,npts);
sigma=1.5; % Adjust this!!

distance=zeros(npts,npts);

for i=1:npts
    i
    for j=1:npts
   % A(i,j)=exp(-norm(STx(i,:)-STx(j,:))^2 / (2*sigma^2)); %% Different distance matrix ? to take into account for single cone response?
    A(i,j) = (STx(i,:)*STx(j,:)')/(norm(STx(i,:)) * norm(STx(j,:)));
    distance(i,j) = norm(STx(i,:)-STx(j,:));
    end
end

figure;
imagesc(A)

figure;
hist(A(:),100)

%A(A<0.93)=0;
A(A<0.5)=0;
A=sparse(A);
figure;
spy(A)

d=sum(A,2);
D=diag(d);
L = eye(npts,npts) - ((inv(D).^(1/2))*A*(inv(D).^(1/2)));
% L =D-A;
L=(L+L')/2;
L=sparse(L);

[E,lam]=eigs(L,20,'sa');
lam=diag(lam);
[lam_sort,idx]=sort(lam,'ascend');

figure;
plot(lam_sort,'*');

E=E(:,idx);

kf=input('Enter number of classes');
U=E(:,1:kf);
for ipt=1:npts
U(ipt,:)=U(ipt,:)/norm(U(ipt,:));
end


[label,C]=kmeans(U,kf);
figure;
plot(U(label==1,1),U(label==1,2),'r*');
hold on
plot(U(label==2,1),U(label==2,2),'b+');
stimdim=size(STx,2);
k_est=zeros(kf,stimdim);

for ik=1:kf
k_est(ik,:)= mean(STx(label==ik,:),1);
end

[label_sort,ii]=sort(label,'ascend');
A_new = A(ii,ii);
figure;
spy(A_new);

figure;
imagesc(A_new);

nc=length(cones);

figure;
icnt=0;
for dim1=1:nc
for dim2=1:nc
     icnt=icnt+1;
    if(dim1~=dim2)
    subplot(nc,nc,icnt);
    plot(X(:,dim1),X(:,dim2),'g.');
    hold on;
    plot(STx(:,dim1),STx(:,dim2),'r.');
    title(sprintf('dim1 %d dim2 %d',dim1,dim2));
    axis square
    end
end
end


su=double((-k_est)==repmat(max(-(k_est)),[kf,1]));
icnt=0;
for isubunits=1:size(su,1)
    if(sum(su(isubunits,:))~=0)
        icnt=icnt+1;
    cells(ref_cell).sub_units(icnt).indices = cones(logical(su(isubunits,:)));
    end
end

       end
       
       
       col='rgbcmykm'
       figure;
       for isu=1:length(cells(ref_cell).sub_units)
           for icone=cells(ref_cell).sub_units(isu).indices
       plot(datarun.cones.centers(icone,1),-datarun.cones.centers(icone,2),strcat(col(isu),'*'));
       axis equal
       hold on
           end
       end

%%
S=double(X'*X);
nc=length(cones);
lam=2;

e=ones(nc,1);
cvx_begin 
variable V(nc,nc) semidefinite
subject to
trace(V)==1
maximize (trace(S*V) - lam*(e'*abs(V)*e))
cvx_end
V

% SparsePCA cvx ? 