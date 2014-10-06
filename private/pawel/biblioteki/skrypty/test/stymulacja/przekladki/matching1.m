sizes=zeros(8,4);

%najpierw n, potem p; najpierw W, potem L

sizes(1,:)=[40 5 200 8];
sizes(2,:)=[40 5 200 8];
sizes(3,:)=[10 5 48 8];
sizes(4,:)=[5 10 18 12];
sizes(5,:)=[2.5 20 6 16];
sizes(6,:)=[1.25 40 2.5 26];
sizes(7,:)=[1 100 1.2 42];
sizes(8,:)=[1 100 1 100];

%sizes(7,:)=[1 100 1.2 42];
%sizes(6,:)=[1.25 40 2 22];
%sizes(5,:)=[2.5 20 6 16];
%sizes(4,:)=[5 10 12 8];
%sizes(3,:)=[10 5 48 8];
%sizes(2,:)=[40 5 150 6];
%sizes(1,:)=[40 5 150 6];

przekl=[1 1 2 1 2 1 2 1];

vn=13.2;
vp=22.7;

volt_shift=zeros(8,3);%dla n, dla p, dla next

%dla zakresu 250u:
volt_shift(2,1)=vn/sqrt(sizes(2,1)*sizes(2,2))*sqrt(2);
volt_shift(2,2)=vp/sqrt(sizes(2,3)*sizes(2,4))*sqrt(2);
volt_shift(2,3)=volt_shift(7,1);

%dla zakresu 1m:
volt_shift(1,1)=vn/sqrt(sizes(1,1)*sizes(1,2))*sqrt(3/2);
volt_shift(1,2)=vp/sqrt(sizes(1,3)*sizes(1,4))*sqrt(3/2);
volt_shift(1,3)=0;

%reszta zakresow:
for i=3:8
    %volt_shift
    v1=vn/sqrt(sizes(i,1)*sizes(i,2));
    v2=vp/sqrt(sizes(i,3)*sizes(i,4));
    if przekl(i)==1
        volt_shift(i,1)=v1*sqrt(3/2);
        volt_shift(i,2)=v2*sqrt(2);
        volt_shift(i,3)=volt_shift(i,1);
    else
        volt_shift(i,1)=v1*sqrt(2);
        volt_shift(i,2)=v2*sqrt(3/2);
        volt_shift(i,3)=volt_shift(i,2);
    end
end

%volt_shift(7,1)=vn/sqrt(sizes(7,1)*sizes(7,2))*sqrt(2);
%volt_shift(7,2)=vp/sqrt(sizes(7,3)*sizes(7,4))*sqrt(2);
%volt_shift(7,3)=volt_shift(7,1);

%volt_shift(8,1)=vn/sqrt(sizes(8,1)*sizes(8,2))*sqrt(3/2);
%volt_shift(8,2)=vp/sqrt(sizes(8,3)*sizes(8,4))*sqrt(3/2);
%volt_shift(8,3)=0;

%cd /home/pawel/pliki/stymulacja/chip/symulacje/2005_01_24;

cd H:\pliki\nauka\stymulacja\chip\symulacje\2005_04_01;

currents={'1m' '250u' '60u' '15u' '4u' '1u' '250n' '60n'};
%res={'1.5k' '6k' '24k' '96k' '380k' '1.5M' '6M' '24M'};
res={'1.45k' '5.8k' '23k' '90k' '360k' '1.44M' '5.8M' '23M'};
l=250;
in=[1:l-1]/l*127;
rozrzut=zeros(8,3,l-1);

for i=1:8    
    cr=currents(i);
    r=res(i);
    
    filename=[cr{1} '_Vg_neg_' r{1} '.out'];
    a=importdata(filename);
    Vg_neg=a(:,2);
    poch_neg=diff(Vg_neg)*1.95*1000; %pochodna wyrazona w milivoltach na LSB
    rozrzut(i,1,:)=volt_shift(i,1)./poch_neg;
    
    filename=[cr{1} '_Vg_pos_' r{1} '.out'];
    a=importdata(filename);
    Vg_pos=a(:,2);
    poch_pos=diff(Vg_pos)*1.95*1000; %pochodna wyrazona w milivoltach na LSB
    rozrzut(i,2,:)=volt_shift(i,2)./poch_pos;
    
    figure(1);
    subplot(2,4,i)
    plot(in,abs(poch_neg'),'bd-',in,abs(poch_pos'),'rd-');
    grid on;

    if przekl(i)==1
        rozrzut(i,3,:)=rozrzut(i,1,:);
    else
        rozrzut(i,3,:)=rozrzut(i,2,:);
    end
end

rozrzut_total=zeros(8,3,l-1);

figure(3)

%dla zakresow od 60u w dol:
rozrzut_input=reshape(rozrzut(2,1,:),1,l-1);
for i=3:8
	%subplot(2,4,i);
	%plot(rozrzut_input);
	%grid on;
    for j=1:l-1
        rozrzut_total(i,1,j)=sqrt(rozrzut_input(j)^2+rozrzut(i,1,j)^2);
	rozrzut_total(i,2,j)=sqrt(rozrzut_input(j)^2+rozrzut(i,2,j)^2);
	if przekl(i)==1
		rozrzut_total(i,3,j)=sqrt(rozrzut_input(j)^2+rozrzut(i,1,j)^2);
		%'n'
	else
		rozrzut_total(i,3,j)=sqrt(rozrzut_input(j)^2+rozrzut(i,2,j)^2);
		%'p'
	end
    end
    subplot(2,4,i);
    a=plot(in,rozrzut_input','b:',in,reshape(rozrzut_total(i,1,:),1,l-1),'r-',in,reshape(rozrzut_total(i,2,:),1,l-1),'k--');
    set(a(1),'LineWidth',1.5);
    set(a(2),'LineWidth',1.5);
    set(a(3),'LineWidth',1.5);
    grid on;
    axis([0 127 0 2]);
    %title=('adfaef')
    xlabel('input DAC code');
    ylabel('sigma of output current spread [LSB]');
    legend('input current','negative current','positive current');

    cr=currents(i)
    %title=('adfaef')
    %filename=[cr{1}% '_Vg_neg_' r{1} '.out'];

    rozrzut_input=reshape(rozrzut_total(i,3,:),1,l-1);    
end
    
%dla 250u oraz 1m:
for j=1:l-1
    rozrzut_total(2,1,j)=rozrzut(2,1,j);
    rozrzut_total(2,2,j)=sqrt(rozrzut(2,1,j)^2+rozrzut(2,2,j)^2);

    rozrzut_total(1,1,j)=rozrzut(1,1,j);
    rozrzut_total(1,2,j)=sqrt(rozrzut(1,1,j)^2+rozrzut(1,2,j)^2);
end

subplot(2,4,1);
a=plot(in,reshape(rozrzut_total(1,1,:),1,l-1),'r-',in,reshape(rozrzut_total(1,2,:),1,l-1),'k--');
set(a(1),'LineWidth',1.5);
set(a(2),'LineWidth',1.5);
grid on;
axis([0 127 0 2]);
xlabel('input DAC code');
ylabel('sigma of output current spread [LSB]');
legend('negative current','positive current');

subplot(2,4,2);
%title('250u');
a=plot(in,reshape(rozrzut_total(2,1,:),1,l-1),'r-',in,reshape(rozrzut_total(2,2,:),1,l-1),'k--');
set(a(1),'LineWidth',1.5);
set(a(2),'LineWidth',1.5);
grid on;
axis([0 127 0 2]);
xlabel('input DAC code');
ylabel('sigma of output current spread [LSB]');
legend('negative current','positive current');
%title('250u');
%for i=5:8
    
