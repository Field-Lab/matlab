function spks=generateResp_LNLN(model,movie,dt,nTrials)

% Get movie at higher resolution
scale = model.gridSzX/size(movie,1);



% weigh inputs for differnet cones
coneInp = zeros(size(movie,3),model.nCones);
coneWindowSz = model.coneWindowSz; gridSzX = model.gridSzX; gridSzY=model.gridSzY;

cone_Map = cell(model.nCones,1);
for icone = 1:model.nCones
coneCenterX = round(model.conesX(icone)); coneCenterY = round(model.conesY(icone));
cone_Map{icone} = model.cone_data{icone}.coneMap(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY));
end

mov_small = cell(model.nCones,1);

for icone=1:model.nCones
mov_small{icone} = zeros(size(cone_Map{1},1),size(cone_Map{1},2),size(movie,3));
end

step=60;
for iframe=1:step:size(movie,3)
    if(rem(iframe,1000)==1)
    iframe
    end
movie_sc = repelem(movie(:,:,iframe:iframe+step-1),scale,scale,1);
for icone = 1:model.nCones   
coneCenterX = round(model.conesX(icone)); coneCenterY = round(model.conesY(icone));
mov_small{icone}(:,:,iframe:iframe+step-1) = movie_sc(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY),:);
end
end


for icone = 1:model.nCones
icone
coneInp(:,icone) = convn(mov_small{icone},repmat(cone_Map{icone}(end:-1:1,end:-1:1),[1,1,1]),'valid');
end

% convolve with tf

movie_c_tf = zeros(size(movie,3),model.nCones);
for icone=1:model.nCones
    cip = coneInp(:,icone);
    ci = zeros(length(cip+30-1),1);
    ci(30:length(cip)+30-1)=cip;
    movie_c_tf(:,icone)=convn(ci,model.ttf,'valid');
end
movie_c_tf = movie_c_tf';

% calculate subunit activation 
SU_inp = model.cone_to_SU_connection*movie_c_tf;
SU_act = model.fs(SU_inp);

% calculate ganglion cell activation
gang_inp = model.SU_gang_weights'*SU_act;
firing_rate = gang_inp * dt;

% generate spikes
spks = poissrnd(repmat(firing_rate,[nTrials,1]));
end