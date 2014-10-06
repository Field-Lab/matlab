cd /home2/pawel/oko/santa_cruz_2003/2003_10_20;
filename='ch28_52uA_01ms_Vpol-11_waterconv';
%!copy filename filename2
%y=convert_data(filename,206,65,20000);

length=60000;

%y=offsets(10,64);

figura=1;
period=20000;
margines=5; 
zakres=2000;  % in miliseconds
stim_channel=28;

style1='bd- ';
style2='go- ';
style3='rx- ';
style4='k*- ';
style5='ch- ';

style6='bd--';
style7='go--';
style8='rx--';
style9='k*--';
style10='ch--';

styles=[style1' style2' style3' style4' style5' style6' style7' style8' style9' style10']';

filename='ch28_11uA_01ms_Vpol-14_waterconv';
y(1,:)=minim(filename,length,stim_channel);

filename='ch28_14uA_01ms_Vpol-14_waterconv';
y(2,:)=minim(filename,length,stim_channel);

filename='ch28_25uA_01ms_Vpol-14_waterconv';
y(3,:)=minim(filename,length,stim_channel);

filename='ch28_52uA_01ms_Vpol-14_waterconv';
y(4,:)=minim(filename,length,stim_channel);

filename='ch28_110uA_01ms_Vpol-14_waterconv';
y(5,:)=minim(filename,length,stim_channel);

filename='ch28_11uA_10ms_Vpol-14_waterconv';
y(6,:)=minim(filename,length,stim_channel);

filename='ch28_14uA_10ms_Vpol-14_waterconv';
y(7,:)=minim(filename,length,stim_channel);

filename='ch28_25uA_10ms_Vpol-14_waterconv';
y(8,:)=minim(filename,length,stim_channel);

filename='ch28_52uA_10ms_Vpol-14_waterconv';
y(9,:)=minim(filename,length,stim_channel);

filename='ch28_110uA_10ms_Vpol-14_waterconv';
y(10,:)=minim(filename,length,stim_channel);

t=[1:64];
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),t,y(5,:),t,y(6,:),t,y(7,:),t,y(8,:),t,y(9,:),t,y(10,:));
legend;

clf
numery=[1 3 5];
subplot(2,1,1);
for i=1:3
	plot(t,y(numery(i),:),styles(i,:));
	styles(i,:)
	hold on;
end
grid on;
axis([1 64 -1200 1200]);
legend('1.1uA, short pulse','2.5uA, short pulse','11uA, short pulse');
ylabel('ADC units');


numery=[6 8 10];
subplot(2,1,2);
for i=1:3
	plot(t,y(numery(i),:),styles(i,:));
	styles(i,:)
	hold on;
end
grid on;
axis([1 64 -1200 1200]);
legend('1.1uA, long pulse','2.5uA, long pulse','11uA, long pulse');
xlabel('channel number');
ylabel('ADC units');
