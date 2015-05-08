%% Analyze spatial STA . 

%% fit GLMs

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
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% in seconds
cellID = 4187;


%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);


%
stim_description=movie_xml;
xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
dashes=find(stim_description=='-');
StimulusPars.type=stim_description(1:dashes(1)-1);
StimulusPars.pixelsize = str2double(stim_description(dashes(1)+1:dashes(2)-1));
StimulusPars.refreshrate = str2double(stim_description(dashes(2)+1:dashes(3)-1));
StimulusPars.RNG = str2double(stim_description(dashes(3)+1:dashes(4)-1));
% try
%     StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
%     StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
% catch
%     StimulusPars.height = 32; 
%     StimulusPars.width = 64;
% end
StimulusPars.tstim = 1/120;
fitframes=stim_length*120; % seconds * 120 frames per second / interval


disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes/StimulusPars.refreshrate);
temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
fitmovie_color=zeros(width,height,3,fitframes);

try
    StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
    StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
catch
    StimulusPars.height =height; 
    StimulusPars.width = width;
end


for i=1:fitframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/StimulusPars.refreshrate));
end

clear temp_fitmovie i 

%% Load spikes 
    master_idx         = find(datarun.cell_ids == cellID);
    
      
        % Spike loading
        spikes=datarun.spikes{master_idx};
    
        % make STA 3D ? 
        %glm_cellinfo.WN_STA = squeeze(sum(glm_cellinfo.WN_STA,3)); % Doubt!!!!!!!
        clear cell_savename
        
        % Align the spikes and the movies;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(datarun.triggers)-1)
            actual_t_start=datarun.triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < datarun.triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
        clear spikes
        spike.home=spikes_adj;
        clear spikes_adj;

        
        %% 
        spksGen = zeros(stim_length*120,1);
        for ispike=1:length(spike.home)
        spksGen(floor(spike.home(ispike)*120)+1)=1;
        end
        mov = fitmovie_color-0.5;
        Filtlen=30;
        
        ifSTC=0;
        [STA,STC] = calculate_sta_stc(spksGen,mov,Filtlen,ifSTC);
        
         [STA_spatial,STC_spatial,u_coll] = calculate_stc_spatial(spksGen,mov,STA,[]);
         
         %%
         
         
         
         k=3;
         [u,s,v]=svd(u_coll_bm(:,1:k)'*u_coll_sm(:,1:k));
         ang_sb=acosd(diag(s))
       
         ang_log=[];
         for irand=1:100
         random_matrix= randn(7,7);
         [u,s,E_rand]=svd(random_matrix);
         [u,s,v]=svd(u_coll_bm(u_coll_bm(:,1)~=0,1:k)'*E_rand(:,1:k));
         ang=acosd(diag(s));
         ang_log(irand)=ang(1);
         end
         
         
         figure;
         hist(ang_log,100);
         hold on;
         plot([ang_sb(1),ang_sb(1)],[0,1],'g');
         
           
         k=3;
         [u,s,v]=svd(u_coll_bbm(:,1:k)'*u_coll_bm(:,1:k));
         ang_bb=acosd(diag(s))
       
         ang_log=[];
         for irand=1:100
         random_matrix= randn(7,7);
         [u,s,E_rand]=svd(random_matrix);
         [u,s,v]=svd(u_coll_bm(u_coll_bm(:,1)~=0,1:k)'*E_rand(:,1:k));
         ang=acosd(diag(s));
         ang_log(irand)=ang(1);
         end
         
         
         figure;
         hist(ang_log,100);
         hold on;
         plot([ang_bb(1),ang_bb(1)],[0,1],'g');
         