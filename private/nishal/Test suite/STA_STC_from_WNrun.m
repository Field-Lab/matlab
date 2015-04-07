function [WNSTA,dimensions,Mask]=STA_STC_from_WNrun(cells, dataset, stim_description, stim_length,dsave,Mask)
% 
% dataset='2014-11-05-2/data009_nps';
% cells={2372,2523}
% stim_description 'RGB-10-2-0.48-11111-32x32'
% stim_length, optional, default 30 min
% dsave, optional, where to save the fit

if nargin==3
    stim_length=1800; 
    d_save='/Volumes/Analysis/nishal/colorglmfits';
elseif nargin==4
    d_save='/Volumes/Analysis/nishal/colorglmfits';
end

%% Stimulus and GLM parameters

% Get stimulus parameters from descriptor
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


% Load datarun
datarun=load_data(dataset);
datarun=load_neurons(datarun);
datarun=load_sta(datarun);
datarun=load_params(datarun);

%% Load Movie
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
clear temp_fitmovie height width i
% testmovie = get_rawmovie(raw_file, testframes);

%% Load Cell Specific Elements   Spikes and STA
for i_cell = 1:length(cells)
    clear glm_cellstruct
    cid = cells{i_cell};
    %[celltype , cell_savename, ~]  = findcelltype(cid, datarun.cell_types);
  
       % Turn RGB movie into greyscale movie
         master_idx         = find(datarun.cell_ids == cid);
         RGB=RGB_weights(datarun,master_idx);
        RGB=RGB/sum(RGB);
        fitmovie=squeeze(RGB(1)*fitmovie_color(:,:,1,:)+ ...
            RGB(2)*fitmovie_color(:,:,2,:)+ ...
            RGB(3)*fitmovie_color(:,:,3,:));
        
        fitmovie=fitmovie-mean(fitmovie(:)); % make it 0 mean!
        
        clear RGB
        
        % Spike loading
        spikes=datarun.spikes{master_idx};
        WNSTA_vision = datarun.stas.stas{master_idx};
        
        
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
        
        %% binary spike train
        binnedResponses = zeros(fitframes,1);
        t_bin        = StimulusPars.tstim;
        home_sptimes = spike.home';
        home_spbins  = ceil(home_sptimes / t_bin);
        binnedResponses(home_spbins)=1;
        
        %% Compute STAs
    movie=fitmovie;
    movieLen=size(movie,3);
    
    movie4D = zeros(size(movie,1),size(movie,2),3,movieLen);
    for iframe=1:size(movie,3);
    movie4D(:,:,1,iframe)=movie(:,:,iframe);
    movie4D(:,:,2,iframe)=movie(:,:,iframe);
    movie4D(:,:,3,iframe)=movie(:,:,iframe);
    end
    mov_params.mov=movie4D;
    
    sta_params.Filtlen=30;

    cell_params.binsPerFrame=1;
    
    response.spksGen=binnedResponses';
    aa=repmat([1:movieLen],[1,1]);
    response.mov_frame_number=aa(:);
    
    
    WNSTA=zeros(size(movie,1),size(movie,2),sta_params.Filtlen);

    sta_params.useTrial=1; % Need to use all 5 trials
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
   
  WNSTA=response.analyse.STA;
        
        %% Compute Mask
           % Clip STA
         cell_p.STAlen=sta_params.Filtlen;
         stas{1}=zeros(size(WNSTA,1),size(WNSTA,2),3,cell_p.STAlen);
         for itime = 1:cell_p.STAlen
         stas{1}(:,:,1,itime)=WNSTA(:,:,itime);
         stas{1}(:,:,2,itime)=WNSTA(:,:,itime);
         stas{1}(:,:,3,itime)=WNSTA(:,:,itime);
         end
         
         [new_stas,totalMaskAccept,CellMasks]=clipSTAs(stas,cell_p);
         if(sum(Mask(:))==0)
         Mask = CellMasks{1};
         end
         figure;
         imagesc(Mask);
         axis image;
         
        %% Compute STCs 
        % Have fitmovie,spike,WNSTA,Mask
tic;
reSTA = WNSTA;

Filtdim1=size(reSTA,1);
Filtdim2=size(reSTA,2);
Filtlen=size(reSTA,3);
reSTA = reSTA(logical(repmat(Mask,[1,1,Filtlen])));

%xxreSTA = reshape(reSTA,[Filtdim1*Filtdim2*Filtlen,1]);
xxreSTA=reSTA;
xxreSTA=xxreSTA/norm(xxreSTA(:));

reSTC=zeros(length(xxreSTA));
indx=[59:1:movieLen];
binnedResponses=binnedResponses';
for iTrial=1

binnedResponsesTrial=binnedResponses(iTrial,:);

framesValid = indx(binnedResponsesTrial(indx)>0);

for iframe=framesValid
    iframe
    a = movie(:,:,iframe:-1:iframe-Filtlen+1);
    a=a(logical(repmat(Mask,[1,1,Filtlen])));
  xx= a;
%  mask!

%  xx=  reshape(mov_new2(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
% Don't mask!

xxnew = (xx-((xx'*xxreSTA)*xxreSTA));
  
reSTC=reSTC+ xxnew*xxnew'*binnedResponsesTrial(iframe);
end
end
reSTC=reSTC / (sum(binnedResponses(:))-1);

[u,s,v]=svds(reSTC,100);
figure;
plot(diag(s),'*')
toc;

    
 %% plot points 
h2=figure('Color','w');

component1 = xxreSTA;
component2 = u(:,1);
subplot(3,2,1);
proj_plot_script
title('STA - STC 1');

component1 = xxreSTA;
component2 = u(:,2);
subplot(3,2,2);
proj_plot_script
title('STA - STC 2');

component1 = u(:,1);
component2 = u(:,2);
subplot(3,2,3);
proj_plot_script
title('STC 1 - STC 2');

component1 = u(:,3);
component2 = u(:,2);
subplot(3,2,4);
proj_plot_script
title('STC 3 - STC 2');

component1 = u(:,1);
component2 = u(:,4);
subplot(3,2,5);
proj_plot_script
title('STC 1 - STC 4');

subplot(3,2,6);
plot(diag(s(1:10,1:10)),'*');
title('STC magnitudes');

%%
figure;
xidx=repmat([1:size(Mask,1)]',[1,size(Mask,2)]);
yidx=repmat(1:size(Mask,2),[size(Mask,1),1]);

x_mask=xidx(logical((Mask)));
y_mask=yidx(logical((Mask)));
%iSTC=1;
%uSq = reshape(u(:,iSTC),[length(x_mask),Filtlen]);
uSq = reshape(reSTA,[length(x_mask),Filtlen]);

STC_comp=zeros(size(WNSTA));
for itime=1:Filtlen
maskedU=uSq(:,itime);
for ivalue=1:length(maskedU)
STC_comp(x_mask(ivalue),y_mask(ivalue),itime)=maskedU(ivalue);
end
end

dimensions.Masked.STA=STC_comp;
dimensions.Masked.STAtimeCourse=squeeze(mean(mean(STC_comp,1),2));
subplot(2,6,6+1);
plot(dimensions.Masked.STAtimeCourse);

[V,I]=max(abs(dimensions.Masked.STAtimeCourse));
subplot(2,6,1);
imagesc(STC_comp(:,:,I)');
colormap gray
axis image
title('STA');

for iSTC=1:5

uSq = reshape(u(:,iSTC),[length(x_mask),Filtlen]);

STC_comp=zeros(size(WNSTA));
for itime=1:Filtlen
maskedU=uSq(:,itime);
for ivalue=1:length(maskedU)
STC_comp(x_mask(ivalue),y_mask(ivalue),itime)=maskedU(ivalue);
end
end

dimensions.Masked.STC=cell(6,1);
dimensions.Masked.timeCourse=cell(6,1);
dimensions.Masked.STC{iSTC}=STC_comp;
dimensions.Masked.timeCourse{iSTC}=squeeze(mean(mean(STC_comp,1),2));
subplot(2,6,6+iSTC+1);
plot(dimensions.Masked.timeCourse{iSTC});

[V,I]=max(abs(dimensions.Masked.timeCourse{iSTC}));
subplot(2,6,iSTC+1);
imagesc(STC_comp(:,:,I)');
colormap gray
axis image
color bar
title(sprintf('STC %d',iSTC));

end
%%
dimensions.xxreSTA=xxreSTA;
dimensions.u=u;
dimensions.Mask=Mask;
dimensions.s=s;    
   
    
end
