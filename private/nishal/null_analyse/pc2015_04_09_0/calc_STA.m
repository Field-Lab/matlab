
function WNSTA=calc_STA(movie,spks)
%% Calculate STA
movieLen=size(movie,4);

    mov_params.mov=movie;
    
    sta_params.Filtlen=30;

    cell_params.binsPerFrame=1;
    
    response.spksGen=zeros(1,movieLen);
    response.spksGen(1,spks)=1;
    aa=repmat([1:movieLen],[1,1]);
    response.mov_frame_number=aa(:);
    
  
    sta_params.useTrial=1; % Need to use all 5 trials
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
    WNSTA = response.analyse.STA;
 
%        figure
%          for itime=1:sta_params.Filtlen
%              itime
%          imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
%          caxis([min(WNSTA(:)),max(WNSTA(:))]);
%          colorbar
%          pause(1/120)
%          end
end