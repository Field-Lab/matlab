function []=offplots(filename);

cd H:/pliki/nauka/stymulacja/chip/testy/VOLTAGE_MODE_OFFSET/

%n1=filename;
n1=filename(1:14);
n2=filename(16:18);
n3=filename(6:7);
n4=filename(16:15);

aname=[n1,'4',n2,'_dc0pos.dat']
a=importdata(aname)';
bname=[n1,'4',n2,'_dc50pos.dat'];
b=importdata(bname)';
cname=[n1,'4',n2,'_dc100pos.dat'];
c=importdata(cname)';
dname=[n1,'4',n2,'_dc200pos.dat'];
d=importdata(dname)';
ename=[n1,'4',n2,'_dc300pos.dat'];
e=importdata(ename)';
fname=[n1,'4',n2,'_dc400pos.dat'];
f=importdata(fname)';
gname=[n1,'4',n2,'_dc500pos.dat'];
g=importdata(gname)';
hname=[n1,'5',n2,'_dc0neg.dat'];
h=importdata(hname)';
iname=[n1,'5',n2,'_dc50neg.dat'];
i=importdata(iname)';
jname=[n1,'5',n2,'_dc100neg.dat'];
j=importdata(jname)';
kname=[n1,'5',n2,'_dc200neg.dat'];
k=importdata(kname)';
lname=[n1,'5',n2,'_dc300neg.dat'];
l=importdata(lname)';
mname=[n1,'5',n2,'_dc400neg.dat'];
m=importdata(mname)';
nname=[n1,'5',n2,'_dc500neg.dat'];
n=importdata(nname)';

aTDSoffset=[(a(1,1)+a(1,258))/2]
bTDSoffset=[(b(1,1)+b(1,258))/2]
cTDSoffset=[(c(1,1)+c(1,258))/2]
dTDSoffset=[(d(1,1))]%+d(1,258))/2];
eTDSoffset=[(e(1,1)+e(1,258))/2]
fTDSoffset=[(f(1,1)+f(1,258))/2]
gTDSoffset=[(g(1,1)+g(1,258))/2]
hTDSoffset=[(h(1,1)+h(1,258))/2]
iTDSoffset=[(i(1,1)+i(1,258))/2]
jTDSoffset=[(j(1,1)+j(1,258))/2]
kTDSoffset=[(k(1,1)+k(1,258))/2]
lTDSoffset=[(l(1,1)+l(1,258))/2]
mTDSoffset=[(m(1,1)+m(1,258))/2]
nTDSoffset=[(n(1,1)+n(1,258))/2]

dac=[127:-1:-127];

a1=[a(1,2:257)-aTDSoffset];
b1=[b(1,2:257)-bTDSoffset+0.05];
c1=[c(1,2:257)-cTDSoffset+0.1];
d1=[d(1,2:257)-dTDSoffset+0.2];
e1=[e(1,2:257)-eTDSoffset+0.3];
f1=[f(1,2:257)-fTDSoffset+0.4];
g1=[g(1,2:257)-gTDSoffset+0.5];
h1=[h(1,2:257)-hTDSoffset];
i1=[i(1,2:257)-iTDSoffset-0.05];
j1=[j(1,2:257)-jTDSoffset-0.1];
k1=[k(1,2:257)-kTDSoffset-0.2];
l1=[l(1,2:257)-lTDSoffset-0.3];
m1=[m(1,2:257)-mTDSoffset-0.4];
n1=[n(1,2:257)-nTDSoffset-0.5];


a1=[a(1,2:257)-aTDSoffset];
b1=[b(1,2:257)-bTDSoffset];
c1=[c(1,2:257)-cTDSoffset];
d1=[d(1,2:257)-dTDSoffset];
e1=[e(1,2:257)-eTDSoffset];
f1=[f(1,2:257)-fTDSoffset];
g1=[g(1,2:257)-gTDSoffset];
h1=[h(1,2:257)-hTDSoffset];
i1=[i(1,2:257)-iTDSoffset];
j1=[j(1,2:257)-jTDSoffset];
k1=[k(1,2:257)-kTDSoffset];
l1=[l(1,2:257)-lTDSoffset];
m1=[m(1,2:257)-mTDSoffset];
n1=[n(1,2:257)-nTDSoffset];



[slope1,lin_err1,dane1]=dac_lin3(dac,a1)
[slope2,lin_err1,dane2]=dac_lin3(dac,b1);
[slope3,lin_err1,dane3]=dac_lin3(dac,c1);
[slope4,lin_err1,dane4]=dac_lin3(dac,d1);
[slope5,lin_err1,dane5]=dac_lin3(dac,e1);
[slope6,lin_err1,dane6]=dac_lin3(dac,f1);
[slope7,lin_err1,dane7]=dac_lin3(dac,g1);
[slope8,lin_err1,dane8]=dac_lin3(dac,h1)
[slope9,lin_err1,dane9]=dac_lin3(dac,i1);
[slope10,lin_err1,dane10]=dac_lin3(dac,j1);
[slope11,lin_err1,dane11]=dac_lin3(dac,k1);
[slope12,lin_err1,dane12]=dac_lin3(dac,l1);
[slope13,lin_err1,dane13]=dac_lin3(dac,m1);
[slope14,lin_err1,dane14]=dac_lin3(dac,n1);


figure(1);
%a=plot(dac,dane1,'bd-',dac,dane3,'bd-',dac,dane4,'bd-',dac,dane5,'bd-',dac,dane6,'bd-',dac,dane7,'bd-',dac,dane10,'bd-',dac,dane11,'bd-',dac,dane12,'bd-',dac,dane13,'bd-',dac,dane14,'bd-')
a=plot(dac,-dane1,'bd-',dac,-dane3,'bd-',dac,-dane4,'bd-',dac,-dane5,'bd-',dac,-dane6,'bd-',dac,-dane7,'bd-',dac,-dane10,'bd-',dac,-dane11,'bd-',dac,-dane12,'bd-',dac,-dane13,'bd-',dac,-dane14,'bd-')
for i=1:11
    set(a(i),'MarkerSize',5);
end
fs=30;
%title1=['volt ',n3,' 120k ',n4,' c1 with offset']; 
%title(title1,'FontSize',18,'FontWeight','demi');
%xlabel('DAC value','FontSize',fs,'FontWeight','demi');
%ylabel('Voltage [V]','FontSize',fs,'FontWeight','demi');
grid on;
%axis([-130 130 -1 1]);
h=gca;
set(h,'FontSize',fs);
set(h,'LineWidth',1.5);
set(h,'XLim',[-140,140],'YLim',[-1 1]);
%axis([-130 130 -1 1]);
clear k;
[P,S]=POLYFIT(dac,dane1,1);
k(1)=P(1,1);
[P,S]=POLYFIT(dac,dane3,1);
k(2)=P(1,1);
[P,S]=POLYFIT(dac,dane4,1);
k(3)=P(1,1);
[P,S]=POLYFIT(dac,dane5,1);
k(4)=P(1,1);
[P,S]=POLYFIT(dac,dane6,1);
k(5)=P(1,1);
[P,S]=POLYFIT(dac,dane7,1);
k(6)=P(1,1);
[P,S]=POLYFIT(dac,dane10,1);
k(7)=P(1,1);
[P,S]=POLYFIT(dac,dane11,1);
k(8)=P(1,1);
[P,S]=POLYFIT(dac,dane12,1);
k(9)=P(1,1);
[P,S]=POLYFIT(dac,dane13,1);
k(10)=P(1,1);
[P,S]=POLYFIT(dac,dane14,1);
k(11)=P(1,1);
mean(k)

figure(2)
%a=plot(-500,-slope1,'bd-',-400,-slope3,'bd-',-300,-slope4,'bd-',-200,-slope5,'bd-',-100,-slope6,'bd-',0,-slope7,'bd-',100,-slope10,'bd-',200,-slope11,'bd-',300,-slope12,'bd-',400,-slope13,'bd-',500,-slope14,'bd-');
plot([-500 -400 -300 -200 -100 0 100 200 300 400 500]', (k-mean(k))./mean(k))
size(k)
figure(4)
a=plot([-500 -400 -300 -200 -100 0 100 200 300 400 500]',-k*1000,'bd-')
set(a,'LineWidth',2);
set(a,'MarkerSize',10);
axis([-550 550 3.08 3.15]);
h=gca;
set(h,'FontSize',30);
set(h','LineWidth',2);
set(h,'XTick',[-400 -200 0 200 400]);
grid on;
figure(3)
dac=[1:256];
a=plot(dac,a1,'bd-',dac,c1,'bd-',dac,d1,'bd-',dac,e1,'bd-',dac,f1,'bd-',dac,f1,'bd-',dac,j1,'bd-',dac,k1,'bd-',dac,l1,'bd-',dac,m1,'bd-',dac,n1,'bd-');

cd pictures;
filename2=[filename,'with_offset'];
print('-dtiff', filename2);
cd ..;