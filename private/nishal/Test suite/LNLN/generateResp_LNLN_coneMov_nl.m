function [spks,SU_inp,movie_c_tf]=generateResp_LNLN_coneMov_nl(model,coneInp,dt,nTrials)

% % Get movie at higher resolution
% scale = model.gridSzX/size(movie,1);
% 
% 
% 
% % weigh inputs for differnet cones
% coneInp = zeros(size(movie,3),model.nCones);
% coneWindowSz = model.coneWindowSz; gridSzX = model.gridSzX; gridSzY=model.gridSzY;
% 
% cone_Map = cell(model.nCones,1);
% for icone = 1:model.nCones
% coneCenterX = round(model.conesX(icone)); coneCenterY = round(model.conesY(icone));
% cone_Map{icone} = model.cone_data{icone}.coneMap(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY));
% 
% end
% 
% mov_small = cell(model.nCones,1);
% 
% step=60;
% for iframe=1:step:size(movie,3)
%     if(rem(iframe,100)==1)
%     iframe
%     end
% movie_sc = repelem(movie(:,:,iframe:iframe+step-1),scale,scale,1);
% for icone = 1:model.nCones   
% coneCenterX = round(model.conesX(icone)); coneCenterY = round(model.conesY(icone));
% mov_small{icone}(:,:,iframe:iframe+step-1) = movie_sc(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY),:);
% end
% end
% 
% for icone = 1:model.nCones
% icone
% mov_small{icone} =gather(mov_small{icone});
% coneInp(:,icone) = convn(mov_small{icone},repmat(cone_Map{icone}(end:-1:1,end:-1:1),[1,1,1]),'valid');
% end

% convolve with tf

movie_c_tf = zeros(size(coneInp,1),model.nCones);
downsamplert = round(dt/(1/120));
newLen = floor(length(model.ttf)/downsamplert);

for icone=1:model.nCones
    cip = coneInp(:,icone);
    ci = zeros(length(cip+newLen-1),1);
    ci(newLen:length(cip)+newLen-1)=cip;
    ttf=downsample(model.ttf,downsamplert);
    movie_c_tf(:,icone)=convn(ci,ttf,'valid');
end
movie_c_tf = movie_c_tf';

% calculate subunit activation 
SU_inp = model.cone_to_SU_connection*movie_c_tf;
SU_act = model.fs(SU_inp);

% calculate ganglion cell activation
gang_inp = model.g(model.SU_gang_weights'*SU_act);
firing_rate = gang_inp * dt;

% generate spikes
spks = poissrnd(repmat(firing_rate,[nTrials,1]));
spks= double(spks~=0); % make 0 or 1 spike

end