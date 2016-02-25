function elecs = getElectrodes_simulation(coneSpacing,arrSz,coneLatticeOrientation)
d = coneSpacing;

xdim = arrSz;
ydim = arrSz;

% Generate lattice vertices
conesX = [];
conesY = [];

% Rotate lattice
theta = coneLatticeOrientation;
R = [cos(theta),-sin(theta);
    sin(theta),cos(theta)];

for ix = 0:xdim
    for iy=-ydim/3:ydim*2/3
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
% elecs.xdim = xdim;
% elecs.ydim = ydim;
end