
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

if(exist('fitmovie_color','var')~=1) % variable fit movie_color does not exist
    
disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,dur,ref] = get_movie(xml_file, datarun.triggers, fitframes/StimulusPars.refreshrate);
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
end
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
        spksGen = spksGen(1:stim_length*120);
        
        spksGen_hr = zeros(stim_length*1200,1);
        for ispike=1:length(spike.home)
            spksGen_hr(floor(spike.home(ispike)*1200)+1)=1;
        end
        spksGen_hr = spksGen_hr(1:stim_length*1200);
         %% find significant mask - 1
        idx=1:length(datarun.cell_ids);
        matlab_cell_ids = idx(datarun.cell_ids==cellID);
        stas=datarun.stas.stas(matlab_cell_ids);
        
        % Load STAs
        
        stas_new=cell(length(stas),1);
        for icell=1:length(stas)
            st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
            for itime=1:size(stas{1},4)
                st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
            end
            stas_new{icell}=st_temp;
        end
        stas=stas_new;
        
        % Used in movie post process
        cell_params.STAlen=30;
        [stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs(stas,cell_params);
        

% Calculate STA
        STA_depth=user_STA_depth;
       
        mov = fitmovie_color-0.5;
        dim1=size(mov,1);
        dim2=size(mov,2);
        iisp = 1:size(mov,4);
        iisp = iisp(spksGen~=0 & iisp' > STA_depth+1);
        
        STA_recalc = zeros(dim1,dim2,3,STA_depth);
        for ilen=1:length(iisp)
        STA_recalc = STA_recalc + mov(:,:,:,iisp(ilen)-(STA_depth-1):iisp(ilen));
        end
        STA_recalc =STA_recalc/numel(iisp);
        
        % Find the mask using freshly calculated STA
        hhh=figure; 
        for itime=1:STA_depth
        subplot(1,2,1);
        imagesc(mean(STA_recalc(:,:,:,itime),3)');
        colormap gray
       caxis([min(STA_recalc(:)),max(STA_recalc(:))]);
       
        axis image
        subplot(1,2,2);
        imagesc(totalMaskAccept2);
      %  colormap gray
        axis image
        %tf 
        pause(0.1); 
        end
            
      ttf=squeeze(mean(mean(mean(STA_recalc.*repmat(totalMaskAccept2,[1,1,3,STA_depth]),3),1),2));
      %tf=tf/max(abs(tf));
      %idx=1:length(tf);
      ttf=1000*ttf(end:-1:1);
      %tf=tf.*double(idx<15)';
     
      figure;
      plot(ttf);
      title('Is this STA depth correct? ');
      pause(0.5);
      
      [maxV,maxSTATime] = max(squeeze(ttf));
      
        [r,c] = find(totalMaskAccept2>0);
        x_coord = [max(min(r)-3,1):min(max(r)+3,size(totalMaskAccept2,1))];
        y_coord = [max(min(c)-3,1):min(max(c)+3,size(totalMaskAccept2,2))];
    
     
 mov=squeeze(mean(mov,3));
 maskedMovdd= filterMov(mov,totalMaskAccept2,squeeze(ttf));
 maskedMovdd2=[maskedMovdd;ones(1,size(maskedMovdd,2))];
 totalMaskAccept = totalMaskAccept2;
 
 clear mov
% [fitGMLM,output] = fitGMLM_afterSTC_simplified(binnedResponses,maskedMov,7,4);
%  
%  x_coord=[26-13:26+13]; y_coord = [26-13:26+13];
% %% Load old GLM fits
% fitGLM = load (cell_glm_fit);
% fitGLM=fitGLM.fittedGLM;
% % tf = fitGLM.linearfilters.Stimulus.time_rk1;
% % figure;
% % plot(tf);
% x_coord = fitGLM.linearfilters.Stimulus.x_coord;
% y_coord = fitGLM.linearfilters.Stimulus.y_coord;
% totalMaskAccept = zeros(40,40);
% totalMaskAccept(x_coord,y_coord)=1;
% 
%  
%         mov = fitmovie_color-0.5;
%         iisp = 1:size(mov,4);
%         iisp = iisp(spksGen~=0 & iisp' > 31);
%         STA_recalc = zeros(40,40,3,30);
%         for ilen=1:length(iisp)
%         STA_recalc = STA_recalc + mov(:,:,:,iisp(ilen)-29:iisp(ilen));
%         end
%         STA_recalc =STA_recalc/numel(iisp);
%         
%         figure; 
%         subplot(1,2,1);
%         imagesc(mean(STA_recalc(:,:,:,26),3));
%         colormap gray
%         subplot(1,2,2);
%         imagesc(totalMaskAccept);
%         colormap gray
%         
%         figure;
%         imagesc(mean(STA_recalc(x_coord,y_coord,:,26),3));colormap gray
%         
% mov=squeeze(mean(mov,3));
% maskedMovdd= filterMov(mov,totalMaskAccept,squeeze(tf));
% maskedMov2dd=[maskedMovdd;ones(1,size(maskedMo