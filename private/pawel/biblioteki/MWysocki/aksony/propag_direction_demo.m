

% jak odczytac dane:
DataPath='scan_new';
PatternNumber=498 
MovieNumber=63
[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
s=reshape(mean(DataTraces0),512,140);
tmp = mean(s(:,135:end),2);
s = s - repmat(tmp,1,140);
st = s(:,7:40);
%plot(st')
%axis([0 40 -300 300])
[minimus,maxtimes] = min(st,[],2);
minimus = -minimus;
%scatter(maxtimes,maxis);

% jak odczytqc wspolrzedne dla konkretnej elektrody
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
%Radius=2;
%Neighbors=electrodeMap.getAdjacentsTo(PatternNumber,Radius)';
%ElectrodeNumber=512;
elec = zeros(512,3); % x,y,maximum spike
for i = 1:512
    elec(i,1)=electrodeMap.getXPosition(i);
    elec(i,2)=electrodeMap.getYPosition(i);
end
elec(:,3) = minimus;
figure(1)
%elec(Neighbors,3) = 1;
elec = [(1:512)' elec];
elec(420,:) = []; % BAD electrode
%scatter(elec(:,1),elec(:,2),elec(:,3),'fill');
%hold on
f = TriScatteredInterp(elec(:,2:3),elec(:,4),'natural');
tic = 10;
xx = -945:tic:945;
yy = -450:tic:450;
[qx,qy] = meshgrid(xx,yy);
%figure
z = f(qx,qy);
surf(qx,qy,z);
%hold off

figure(2)

maksimum = max(max(z));
dist = 150;
[vert_max,xpos] = max(z,[],2);
[hor_max,ypos] = max(z,[],1);
% pionowo
xvert = xx(xpos);
yvert = yy;
yvert(vert_max < maksimum/7) = [];
xvert(vert_max < maksimum/7) = [];
check = checkDistance(xvert,yvert,elec(elec(:,1)==PatternNumber,2),elec(elec(:,1)==PatternNumber,3),dist);
yvert(check) = [];
xvert(check) = [];
scatter(xvert,yvert,'y');
hold on
xhor = xx;
yhor = yy(ypos);
xhor(hor_max < maksimum/7) = [];
yhor(hor_max < maksimum/7) = [];
check = checkDistance(xhor,yhor,elec(elec(:,1)==PatternNumber,2),elec(elec(:,1)==PatternNumber,3),dist);
xhor(check) = [];
yhor(check) = [];
scatter(xhor,yhor,'r');
axis([-1000 1000 -500 500]);
hold off
pion = corrcoef(xvert,yvert);
poziom = corrcoef(xhor,yhor);
if (pion(1,2) < poziom(1,2))
    points = [xvert' yvert'];
    orient = 1;
else
    points = [xhor' yhor'];
    orient = 0;
end
figure(10)
subplot(2,2,1)
scatter(points(:,1),points(:,2),'b');
hold on
axis([-1000 1000 -500 500]);
wsp = polyfit(points(:,1),points(:,2),1);

edge_points = zeros(4,3);
% right edge
x = 945;
y = wsp(1)*x + wsp(2);
edge_points(1,1) = x;
edge_points(1,2) = y;
if y >= -450 && y <= 450
    edge_points(1,3) = 1;
end
% bottom edge
y = -450;
x = (y - wsp(2))/wsp(1);
edge_points(2,1) = x;
edge_points(2,2) = y;
if x >= -945 && x <= 945
    edge_points(2,3) = 1;
end
% left edge
x = -945;
y = wsp(1)*x + wsp(2);
edge_points(3,1) = x;
edge_points(3,2) = y;
if y >= -450 && y <= 450
    edge_points(3,3) = 1;
end
% top edge
y = 450;
x = (y - wsp(2))/wsp(1);
edge_points(4,1) = x;
edge_points(4,2) = y;
if x >= -945 && x <= 945
    edge_points(4,3) = 1;
end
edge_points(edge_points(:,3)==0,:) = [];
edge_points(:,3) = [];
plot(edge_points(:,1),edge_points(:,2),'r');
hold off

for fig = 2:4
    v = [edge_points(2,1)-edge_points(1,1);edge_points(2,2)-edge_points(1,2)];
    v_norm = sqrt(v'*v)
    v = v/v_norm;
    obrot = [0 -1; 1 0];
    vp = obrot*v;

    rmax = 200;
    dr = 10;
    v0 = edge_points(1,:)';
    points = zeros(floor(v_norm/tic)+1,3);
    for d = 0:1:floor(v_norm/tic)
       sect = getIntersection( v0 + 10*d*v, vp, rmax, dr, f );
       points(d+1,:) = sect(sect(:,4)==max(sect(:,4)),2:4);
    end
    points(points(:,3) < maksimum/7,:) = [];
    check = checkDistance(points(:,1),points(:,2),elec(elec(:,1)==PatternNumber,2),elec(elec(:,1)==PatternNumber,3),dist);
    points(check,:) = [];
    subplot(2,2,fig)
    scatter(points(:,1),points(:,2),'b');
    hold on
    axis([-1000 1000 -500 500]);

    wsp = polyfit(points(:,1),points(:,2),1);

    edge_points = zeros(4,3);
    % right edge
    x = 945;
    y = wsp(1)*x + wsp(2);
    edge_points(1,1) = x;
    edge_points(1,2) = y;
    if y >= -450 && y <= 450
        edge_points(1,3) = 1;
    end
    % bottom edge
    y = -450;
    x = (y - wsp(2))/wsp(1);
    edge_points(2,1) = x;
    edge_points(2,2) = y;
    if x >= -945 && x <= 945
        edge_points(2,3) = 1;
    end
    % left edge
    x = -945;
    y = wsp(1)*x + wsp(2);
    edge_points(3,1) = x;
    edge_points(3,2) = y;
    if y >= -450 && y <= 450
        edge_points(3,3) = 1;
    end
    % top edge
    y = 450;
    x = (y - wsp(2))/wsp(1);
    edge_points(4,1) = x;
    edge_points(4,2) = y;
    if x >= -945 && x <= 945
        edge_points(4,3) = 1;
    end
    edge_points(edge_points(:,3)==0,:) = [];
    edge_points(:,3) = [];
    plot(edge_points(:,1),edge_points(:,2),'r');
    hold off
end

odl = [160 200 300 400 500 600];
r = -rmax:dr:rmax;
wykresy = zeros(length(odl),length(r),2);
przes = [-wsp(2)/wsp(1);0]; % x = -b/a, y = 0

for i = 1:length(odl)
    v0 = v*(v'*(elec(elec(:,1)==PatternNumber,2:3)' - przes)) + odl(i)*v + przes;
    for j = 1:length(r)
       vr = v0 + r(j)*vp;
       wykresy(i,j,1) = r(j);
       wykresy(i,j,2) = f(vr(1),vr(2));
    end
end

figure(4)
subplot(2,1,1);
scatter(elec(:,2),elec(:,3),elec(:,4),'fill');
axis equal
axis([-1000 1000 -500 500]);
hold on
%scatter(elec(:,1),elec(:,2),elec(:,3),'fill');
%hold on
scatter(v0(1),v0(2),'r','fill');
inters(:,1) = v0 + rmax*vp;
inters(:,2) = v0 - rmax*vp;
plot(inters(1,:),inters(2,:),'r');
plot(edge_points(:,1),edge_points(:,2),'r');
scatter(elec(elec(:,1)==PatternNumber,2),elec(elec(:,1)==PatternNumber,3),'y','fill');
hold off
subplot(2,1,2);
plot(wykresy(end,:,1),wykresy(end,:,2));