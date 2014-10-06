filename2=[filename 'conv']
%!copy filename filename2
%y=convert_data(filename,206,65,20000);

length=60000;

figura=1;
period=20000;
margines=5; 
zakres=2000;  % in miliseconds
y=stim_analyze_bez_usr3(filename2,length,period,stim_channels,figura,margines,zakres,1)


figura=2;
period=20000;
margines=5; 
zakres=40;  % in miliseconds
y=stim_analyze_bez_usr3(filename2,length,period,stim_channels,figura,margines,zakres,1)
