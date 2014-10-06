cd I:\dane\2004\26maja\2003-12-08-0-001
%cd /mnt/data4/dane/2004/26maja/2003-12-08-0-001;
filenames={'electrode45.txt' 'electrode137.txt' 'electrode148.txt' 'electrode155.txt' 'electrode169.txt' 'electrode221.txt' 'electrode298.txt' 'electrode362.txt' 'electrode416.txt' 'electrode499.txt'};

%parametry detekcji - b.wazne, aby nie zmieniac!!! bo wtedy numery
%znalezionych spikow sie nie zgodza; do szukanai spikow: funkcja
%oversamp_report_onechannel_guinepig
prog=60; 
histereza=30;
znak=-1;
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);

N=2;
filtr_dl=20;
freq_gr=0.98;
filter_param=struct('N',N,'order',filtr_dl,'freq',freq_gr);

margins=[25 50];

[y1,y2]=inveyefilter_hayes(200,50,50,4000,20000);
RC_filter=y2;

figures=[1 2 3];
%1. Spiki typu "cell body":
%I:\dane\2004\26maja\2003-12-08-0-001 - kolejne pliki:
%electrode45: cell body: 85, axonal: 185,300
%electrode137: cell body: 30, 240 axonal:231, 234
%electrode148: cell body 227, axonal: 106, 254
%electrode155: 2,63 axonal: 285
%electrode169: cb: 158, 181 axonal: 5,285
%el221: cb - 10,147 axonal:none
%el298: cb - 1,20; axonal:227, 228
%el362: cb - 4,5; axonal: 3, 225, 227, 281 (mnooostwo!)
%el416: cb - 61,68,167; axonal: 80, 161, 180 (duzo)
%el499: cb - 14, 17; axonal:none
cb_spikes=zeros(19,2);
cb_spikes(:,1)=[1 2 2 3 4 4 5 5 6 6 7 7 8 8 9 9 9 10 10]';
%cb_spikes(:,2)=[85 30 240 227 2 63 158 181 10 147 1 20 4 5 61 68 167 14 17]';
cb_spikes(:,2)=[1132 30 240 227 2 63 158 181 10 147 1 20 4 5 61 68 167 14 17]';

spikes=cb_spikes(1:6,:);
figura=4;
%y=few_spikes_stat(filenames,spikes,detect_param,filter_param,figura,margins);
%break;

y=oversamp_rep_multichns(filenames,spikes,detect_param,filter_param,margins,RC_filter,figures);
ax_spikes(:,1)=[1 1 2 2 3 3 4 5 5 7 7 8 8 8 8 9 9 9]';
ax_spikes(:,2)=[185 300 231 234 106 254 285 5 285 227 228 3 225 227 281 80 161 180]';
figura=5;
spikes=ax_spikes([1 2 5 11 12 14],:);
y=few_spikes_stat(filenames,spikes,detect_param,filter_param,figura,margins);
