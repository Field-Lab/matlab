function []=cgplots(filename);

%H:/pliki/nauka/stymulacja/chip/testy/2006_04_21/
%dane0=-importdata('curr_r0_1200k_b4_c1')';
cd H:/pliki/nauka/stymulacja/chip/testy/CURRENT_MODE/

filename1=[filename,'.dat'];
dane0=importdata(filename1)';

l=length(filename);

rc=filename(1,1:4);
range=filename(1,7);

if(l==20)
    res=filename(1,9:14);
    %filename2=filename(1,1:25);
    bond=filename(1,17);
    chip=filename(1,20);
    if(res=='19200k')
    R=19200e3;
    end
    
end
if(l==19)
    res=filename(1,9:13);
    %filename2=filename(1,1:26);
    bond=filename(1,16);
    chip=filename(1,19);
    if(res=='1200k')
    R=1200e3;
    end
    if(res=='4700k')
    R=4700e3;
    end
end
if(l==18)
    res=filename(1,9:12);
    %filename2=filename(1,1:27);
    bond=filename(1,15)
    chip=filename(1,18);
    if(res=='1190')
    R=1190;
    end
    if(res=='4800')
    R=4800;
    end
    if(res=='299k')
    R=299e3;
    end
end
if(l==17)
    res=filename(1,9:11);
    bond=filename(1,14);
    chip=filename(1,17);
    if(res=='306') 
    R=306;
    end
    if(res=='19k')
    R=19.1e3;
    end
    if(res=='75k')
    R=75e3;
    end
end
    
dac=[-127:127];

%dane1=[dane0(1,2:257)];

%TDSoffset=[(dane0(1,1)+dane0(1,258))/2];

[slope1,lin_err1,dane]=dac_lin3(dac,dane0);
slope2=slope1/R;
if(slope2<1e-8)
    slope3=1e9*slope2
    slope4 = mat2str(slope3);
    slope = slope4(1,1:5);
    errtitle=[rc,' r',range,' ',res,' b',bond,' c',chip,' slope=',slope,'nA/LSB'];
end
if(slope2>1e-8)
    slope3=1e6*slope2
    slope4 = mat2str(slope3);
    slope = slope4(1,1:5);
    errtitle=[rc,' r',range,' ',res,' b',bond,' c',chip,' slope=',slope,'uA/LSB'];
end

if(rc=='curr')
    datatitle=[rc,' r',range,' ',res,' b',bond,' chip',chip,' data figure'];
    %errtitle=[rc,' r',range,' ',res,' b',bond,' c',chip,' slope=',slope,'nA/LSB'];
    %else datatitle=['current data r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
    % errtitle=['current linearity error r',range,'res=',res,' b',bond,' chip',chip,' ',offset];
end

figure(1);
plot(dac,dane,'db-');
title(datatitle,'FontSize',22,'FontWeight','demi');
xlabel('DAC (number)','FontSize',22,'FontWeight','demi');
ylabel('Voltage [V]','FontSize',22,'FontWeight','demi');
grid on;
h=gca;
set(h,'FontSize',16);
set(h,'LineWidth',1.5);

cd pictures;
filename2=[filename,'_data'];
print('-dtiff', filename2);
cd ..;

%dane2=dane1-TDSoffset;
[slope1,lin_err1,dane]=dac_lin3(dac,dane0);
%slope1;

figure(2);
plot(dac,lin_err1,'bd-');
title(errtitle,'FontSize', 22,'FontWeight', 'demi');
xlabel('DAC (number)','FontSize', 22,'FontWeight', 'demi');
ylabel('linearity error [LSB]','FontSize', 22,'FontWeight', 'demi');
h=gca;
set(h,'FontSize',16);
set(h,'LineWidth',1.5);
grid on;
h=gca;

cd pictures;
filename2=[filename,'_linerr'];
print('-dtiff',filename2);
cd ..;
