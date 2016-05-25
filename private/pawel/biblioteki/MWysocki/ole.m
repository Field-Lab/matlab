%close all

% jak odczytac dane:
DataPath='G:\MWysocki\scan_new';
PatternNumber=498 
MovieNumber=53
[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
s=reshape(mean(DataTraces0),512,140);
tmp = mean(s(:,135:end),2);
s = s - repmat(tmp,1,140);
st = s(:,7:40);
figure(10)
plot(st')
axis([0 40 -300 300])
[minimus,maxtimes] = min(st,[],2);
minimus = -minimus;
%scatter(maxtimes,maxis);

% jak odczytqc wspolrzedne dla konkretnej elektrody
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
Radius=2;
Neighbors=electrodeMap.getAdjacentsTo(PatternNumber,Radius)';
%ElectrodeNumber=512;
elec = zeros(512,3); % x,y,maximum spike
for i = 1:512
    elec(i,1)=electrodeMap.getXPosition(i);
    elec(i,2)=electrodeMap.getYPosition(i);
end
elec(:,3) = minimus;
figure(1)
%elec(Neighbors,3) = 1;
scatter(elec(:,1),elec(:,2),elec(:,3),'fill');
hold on
f = TriScatteredInterp(elec(:,1:2),elec(:,3),'natural');
xx = -945:10:945;
yy = -450:10:450;
[qx,qy] = meshgrid(xx,yy);
%figure
z = f(qx,qy);
figure(7)
surf(qx,qy,z);
figure
elec_trun = elec;
elec_trun(Neighbors,3) = NaN;
f2 = TriScatteredInterp(elec_trun(:,1:2),elec_trun(:,3),'natural');

[ymax,xpos] = max(z,[],2);
[xmax,ypos] = max(z,[],1);
%figure
scatter(xx(xpos),yy,'y');
%hold on
scatter(xx,yy(ypos),'r');
hold off
pion = corrcoef(xx(xpos),yy);
poziom = corrcoef(xx,yy(ypos));
if (pion(1,2) > poziom(1,2))
    points = [xx(xpos)' yy'];
else
    points = [xx' yy(ypos)'];
end

figure(2)
surf(qx,qy,f2(qx,qy));

figure(3)
scatter(points(:,1),points(:,2),'k');
axis([-1000 1000 -500 500]);

% yy = -450:60:450;
% points = zeros(length(yy),2);
% for i = 1:length(yy)
%    el = elec(elec(:,2) == yy(i),:);
%    points(i,1) = sum(el(:,1).*el(:,3))/sum(el(:,3));
% end
% points(:,2) = yy';
% figure
% scatter(points(:,1),points(:,2),'r');
% axis([-1000 1000 -500 500]);
% hold on

wsp = polyfit(points(:,1),points(:,2),1);
figure(15)
subplot(2,1,1);
scatter(elec_trun(:,1),elec_trun(:,2),elec_trun(:,3)*0.6,'fill');
axis equal
hold on
y = wsp(1)*xx + wsp(2);
plot(xx,y,'r');
axis([-1000 1000 -500 500]);
v = [xx(end)-xx(1);y(end)-y(1)];
v = v/sqrt(v'*v);
obrot = [0 -1; 1 0];
vp = obrot*v;
% v = 300*v;
% vp = 300*vp;
% figure
% scatter([0 v(1) vp(1)],[0 v(2) vp(2)]);
% axis equal
% axis([-1000 1000 -500 500]);

rmax = 200;
dr = 10;
odl = 780;
przes = [-wsp(2)/wsp(1);0]; % x = -b/a, y = 0
v0 = v*(v'*(elec(PatternNumber,1:2)' - przes)) + odl*v + przes;
r = -rmax:dr:rmax;
wykres = zeros(length(r),2);
for i = 1:length(r)
   vr = v0 + r(i)*vp;
   wykres(i,1) = r(i);
   wykres(i,2) = f2(vr(1),vr(2));
end

%scatter(elec(:,1),elec(:,2),elec(:,3),'fill');
%hold on
scatter(v0(1),v0(2),'r','fill');
inters(:,1) = v0 + rmax*vp;
inters(:,2) = v0 - rmax*vp;
plot(inters(1,:),inters(2,:),'r');
scatter(elec(PatternNumber,1),elec(PatternNumber,2),'y','fill');
hold off
subplot(2,1,2);
plot(wykres(:,1),wykres(:,2),'bd-');