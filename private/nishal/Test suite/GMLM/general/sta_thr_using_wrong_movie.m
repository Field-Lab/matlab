% using wrong movie, get thresold
function threshold_wrongSTA=sta_thr_using_wrong_movie(wrong_xml,datarun,stim_length,cellID)
stim_description=wrong_xml;
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

%if(exist('fitmovie_color_wrong','var')~=1) % variable fit movie_color does not exist
    
disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes/StimulusPars.refreshrate);
temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
fitmovie_color_wrong=zeros(width,height,3,fitframes);

try
    StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
    StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
catch
    StimulusPars.height =height; 
    StimulusPars.width = width;
end


for i=1:fitframes
    fitmovie_color_wrong(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/StimulusPars.refreshrate));
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
        spksGen = spksGen(1:stim_length*120);
        
        spksGen_hr = zeros(stim_length*1200,1);
        for ispike=1:length(spike.home)
            spksGen_hr(floor(spike.home(ispike)*1200)+1)=1;
        end
        spksGen_hr = spksGen_hr(1:stim_length*1200);
         %% calculate STA
             STA_depth=60;
       
        
        mov = fitmovie_color_wrong-0.5;
         dim1=size(mov,1);
        dim2=size(mov,2);
        iisp = 1:size(mov,4);
        iisp = iisp(spksGen~=0 & iisp' > STA_depth+1);
        
        STA_recalc = zeros(dim1,dim2,3,STA_depth);
        for ilen=1:length(iisp)
        STA_recalc = STA_recalc + mov(:,:,:,iisp(ilen)-(STA_depth-1):iisp(ilen));
        end
        STA_recalc =STA_recalc/numel(iisp);
     
        
        
        %% 
        threshold_wrongSTA = sig_stixels_threshold(STA_recalc, 0.6);
        
end