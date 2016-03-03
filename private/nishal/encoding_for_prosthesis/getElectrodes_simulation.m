function elecs = getElectrodes_simulation(coneSpacing,arrSz,coneLatticeOrientation,gridSz)

gridSzX = gridSz;
gridSzY = gridSz;

d = coneSpacing;
coneGaussSd = coneSpacing/4;

xdim = arrSz;
ydim = arrSz;

% Generate lattice vertices
conesX = [];
conesY = [];

% Rotate lattice
theta = coneLatticeOrientation;
R = [cos(theta),-sin(theta);
    sin(theta),cos(theta)];

for ix = round(-0.5*xdim):xdim
    for iy=-ydim/3:ydim*3/3
        conx = (ix-1)*sqrt(3)*d/2+10;
        cony = (d/2)*rem(ix,2) + (iy-1)*d;
        % rotate
        rot = R*[conx;cony];
        %rot = rot+randn(2,1)*coneLocSd;
        
        % jitter
        conx=rot(1);cony=rot(2);
        
        %if(conx <=gridSzX & cony <=gridSzY & conx>=1 & cony>=1)
        conesX = [conesX;conx];
        conesY = [conesY;cony];
        %end
        
    end
end

nCones = size(conesX,1);

elecs.x = conesX;
elecs.y = conesY;
elecs.nCones=nCones;

% elecs.xdim = xdim;
% elecs.ydim = ydim;

%% 
suX = [];
suY = [];
for ix = round(-0.5*xdim):xdim
    for iy=-ydim/3:ydim*3/3
        conx = (ix-1)*sqrt(3)*d/2 + sqrt(3)*d/(4) +10;
        cony = (iy-1)*d/2;
        % rotate
        rot = R*[conx;cony];
        %rot = rot+randn(2,1)*coneLocSd;
        
        % jitter
        conx=rot(1);cony=rot(2);
        
        %if(conx <=gridSzX & cony <=gridSzY & conx>=1 & cony>=1)
        suX = [suX;conx];
        suY = [suY;cony];
        %end
        
    end
end

elecs.suX=suX;
elecs.suY=suY;

elecs.nSU = length(elecs.suX);

%% assign electrodes to electrode SU
su_elec = sparse(elecs.nSU,elecs.nCones);
for isu=1:elecs.nSU
    if(rem(isu,100)==1) 
        isu 
    end
dists = (elecs.x-elecs.suX(isu)).^2 + (elecs.y-elecs.suY(isu)).^2;
[v,idx] =sort(dists,'ascend');
su_elec(isu,idx(1:3))=rand(1,3);
end
elecs.su_elec = su_elec;
% 
% figure;
% plot(elecs.x,elecs.y,'r.');
% hold on;
% plot(elecs.suX,elecs.suY,'b.');
% hold on;
% for isu=1:elecs.nSU
% iidx = su_elec(isu,:);
% xx = elecs.x(logical(iidx));
% yy = elecs.y(logical(iidx));
% aa = convhull([xx,yy]);
% plot(xx(aa),yy(aa));
% hold on;
% end

%% 
% Make cone input weights
cone_data=cell(nCones,1);
totalConeMap = sparse(gridSzX,gridSzY);

stas = zeros(nCones,numel(totalConeMap));

for icone = 1:nCones
coneMap = zeros(gridSzX,gridSzY);
coneCenterX = round(conesX(icone));
coneCenterY = round(conesY(icone));
if(coneCenterX>=1 & coneCenterX<=gridSzX & coneCenterY>=1 & coneCenterY<=gridSzY )
coneMap(coneCenterX,coneCenterY) = 1;
coneMap = imgaussfilt(coneMap , coneGaussSd^2);
coneMap=coneMap/max(coneMap(:));

coneMap=coneMap(1:gridSzX,1:gridSzY);

%coneWindow = 0*coneMap;
%coneWindow(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY))=1;
end
%cone_data{icone}.coneWindow = logical(coneWindow);
cone_data{icone}.coneMap = coneMap;

stas(icone,:) = coneMap(:);
totalConeMap = totalConeMap + coneMap;
end


figure;
imagesc(totalConeMap');
colormap gray
title('Cone map');

%% 

elecs.gridSzX = gridSzX;
elecs.gridSzY = gridSzY;
elecs.coneSpacing=coneSpacing;
elecs.coneGaussSd=coneGaussSd;
elecs.cone_data=cone_data;
elecs.totalConeMap = totalConeMap;
elecs.stas=stas;
end