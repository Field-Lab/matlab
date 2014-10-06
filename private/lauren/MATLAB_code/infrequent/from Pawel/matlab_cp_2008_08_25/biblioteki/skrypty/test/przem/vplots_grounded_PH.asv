function []=vgplots(filename);

cd H:/pliki/nauka/stymulacja/chip/testy/VOLTAGE_MODE_GROUND

filename1=[filename,'.dat'];

dane0=importdata(filename1)';
dane1=fliplr(dane0);                    %funkcja fliplr odwraca macierz ze strony prawej na lewa

l=length(filename);

vc=filename(1,1:4);
range=filename(1,7);
res=filename(1,9:12);
bond=filename(1,15);
chip=filename(1,18);
if(l==29)
    offset=filename(1,20:25);
    filename2=filename(1,1:25);
end
if(l==25)
    offset=filename(1,20:25);
    filename2=filename(1,1:25);
end
if(l==31)
    offset=filename(1,20:27);
    filename2=filename(1,1:27);
end

dac=[-127:127];

%dane1=[dane0(1,1:256)];

%TDSoffset=[(dane0(1,1)+dane0(1,258))/2];

[slope1,lin_err1,dane]=dac_lin4(dac,dane0);
slope1
%slope2 = 1000*slope1;
%slope3 = mat2str(slope2);
%slope=slope3(1,1:6);
phisoff=1000*dane(1,128);
phisoff2=mat2str(phisoff);
phisoff3=phisoff2(1,1:5);
if(vc=='volt')
   % datatitle=[vc,' r',range,' ',res,' b',bond,' chip',chip,' offset=',phisoff3,'mV data figure'];
    %errtitle=[vc,' r',range,' ',res,' b',bond,' chip',chip,' grd. slope',slope,'mV/LSB'];
    %else datatitle=['current data r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
    % errtitle=['current linearity error r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
end

figure(1);
datatitle=['   ',vc,' r',range,' ',res,' b',bond,' chip',chip,' offset=',phisoff3,'mV data figure'];
plot(dac,dane,'db-');
%title(datatitle,'FontSize',22,'FontWeight','demi');
xlabel('DAC setting','FontSize',22,'FontWeight','demi');
ylabel('output voltage [V]','FontSize',22,'FontWeight','demi');
grid on;
axis([-140 140 -2 2]);
h=gca;
set(h,'FontSize',20);
set(h,'LineWidth',2);

cd pictures;
filename3=[filename2,'_data'];
print('-dtiff', filename3);
cd ..;

%dane2=dane1-TDSoffset;
[slope1,lin_err1,dane]=dac_lin4(dac,dane1);
slope2 = 1000*slope1;
slope3 = mat2str(slope2);
slope=slope3(1,1:6);
slope1
slope3

errtitle=[vc,' r',range,' ',res,' b',bond,' chip',chip,' grd. slope',slope,'mV/LSB'];
figure(2);
as=plot(dac,lin_err1,'bd-');
%set(as,'MarkerSize',10);
%title(errtitle,'FontSize', 22,'FontWeight', 'demi');
%xlabel('DAC (number)','FontSize', 20,'FontWeight', 'demi');
%ylabel('linearity error [LSB]','FontSize', 20,'FontWeight', 'demi');
h=gca;
set(h,'FontSize',30);
set(h,'LineWidth',2);
axis([-140 140 -0.5 1.5]);
grid on;

cd pictures;
filename3=[filename2,'_linerr'];
print('-dtiff',filename3);
cd ..;

figure(11)
plot(dac,dane)