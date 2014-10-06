function []=voffplots_dokt(filename);

cd H:/pliki/nauka/stymulacja/chip/testy/VOLTAGE_MODE_OFFSET
%H:/pliki/nauka/stymulacja/chip/testy/2006_04_21/
filename1=[filename,'.dat'];
dane0=importdata(filename1)';
%dane0=-importdata('volt_r3_120k_b5_c1_dc0neg.dat')';

l=length(filename);

vc=filename(1,1:4);
range=filename(1,7);
res=filename(1,9:12);
bond=filename(1,15);
chip=filename(1,18);
if(l==25)
    offset=filename(1,20:25);
    filename2=filename(1,1:25);
end
if(l==26)
    offset=filename(1,20:26);
    filename2=filename(1,1:26);
end
if(l==27)
    offset=filename(1,20:27);
    filename2=filename(1,1:27);
end

dac=[-127:127];

dane1=[dane0(1,2:257)];

TDSoffset=[(dane0(1,1)+dane0(1,258))/2];

[slope1,lin_err1,dane]=dac_lin3(dac,dane1);
slope = mat2str(slope1);                        %konwersja macierzy na string

if(vc=='volt')
    datatitle=[vc,' r',range,' ',res,' b',bond,' chip',chip,' ',offset,' data figure'];
    errtitle=[vc,' r',range,' ',res,' b',bond,' chip',chip,' ',offset,' ',slope,'aaa'];
else datatitle=['current data r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
     errtitle=['current linearity error r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
end

figure(1);
plot(dac,dane,'db-');
title(datatitle,'FontSize',18,'FontWeight','demi');
xlabel('DAC (number)','FontSize',18,'FontWeight','demi');
ylabel('Voltage [V]','FontSize',18,'FontWeight','demi');
grid on;
h=gca;
set(h,'FontSize',16);

cd pictures;
filename3=[filename2,'_data'];
print('-dtiff', filename3);
cd ..;

dane2=dane1-TDSoffset;
[slope1,lin_err1,dane]=dac_lin3(dac,dane2);
%slope1;

figure(2);
plot(dac,lin_err1,'bd-');
title(errtitle,'FontSize', 18,'FontWeight', 'demi');
xlabel('DAC (number)','FontSize', 18,'FontWeight', 'demi');
ylabel('linearity error [LSB]','FontSize', 18,'FontWeight', 'demi');
h=gca;
set(h,'FontSize',16);
grid on;

cd pictures;
filename3=[filename2,'_linerr'];
print('-dtiff',filename3);
cd ..;
