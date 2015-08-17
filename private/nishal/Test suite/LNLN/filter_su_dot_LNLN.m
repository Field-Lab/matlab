function [SU_inp]=filter_su_dot_LNLN(model,movie)

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

step=1;
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
    movie_c_tf(:,icone)=cip;
end
movie_c_tf = movie_c_tf';

% calculate subunit activation 
SU_inp = model.cone_to_SU_connection*movie_c_tf;
SU_act = model.fs(SU_inp);
end