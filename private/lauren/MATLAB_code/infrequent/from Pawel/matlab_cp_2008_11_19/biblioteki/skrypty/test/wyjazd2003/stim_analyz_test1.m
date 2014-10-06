%y=convert_data(filename,206,65,20000);
cd C:\stim_data\data\2003_10_17;
filename='ch17_37uA_01ms_water';
%y=convert_data(filename,206,65,20000);
stim_channels=17;

filename=[filename 'conv']
length=800000;

figura=1;
period=4000;
margines=200; 
zakres=200;  % in miliseconds
filename
%y=stim_analyze_bez_usr(filename,length,period,stim_channels,figura,margines,zakres);

figura=2;
period=2000;
margines=10; 
zakres=50;  % in miliseconds
%y=stim_analyze_bez_usr(filename,length,period,stim_channels,figura,margines,zakres,0);

figura=3;
period=40000;
margines=5; 
zakres=10000;  % in miliseconds
y=stim_analyze_bez_usr(filename,length,period,stim_channels,figura,margines,zakres,0);