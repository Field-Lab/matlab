%% 20120902_1
clear
date='20120902_1'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 -1 -1 1 -1 1 -1 -1 1 0 1 -1 -1 1 -1 1 1 0 1 -1 1 -1 0 0 -1 0 1 -1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20120902_2
clear
date='20120902_2'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 1 -1 0 0 0 1 1 1 -1 -1 -1 0 0 1 1 -1 1 0 -1 1 -1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20121023_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 -1 0 -1 -1 0 0 1 0 0 -1 1 0 0 0 1 0 0 1 -1 1 -1 -1 1 1 1 -1 1 -1 0 1 0 1 1 1 1 1 -1 0 1 1 1 -1 1 0 0 0];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20121026_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121026_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[0 -1 1 -1 -1 1 -1 0 1 0 -1 -1 0 1 -1 1 1 -1 1 0 1 -1 -1 1 1 1 0 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20121023
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 1 1 1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20121002
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121002'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
onOff=[-1 1 -1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 -1 -1 -1 1 1 1 -1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20120928
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120928'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
onOff=[-1 -1 1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 1 -1 -1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130301
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 1 -1 1 1 -1 -1 -1 1 1 1 1 1 0 0 1 1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130301_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_1'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 -1 0 1 -1 1 0 0 1 -1 0 -1 0 -1 1 1 0 -1 -1 -1 1 -1 1 1 -1 1 0 1 -1 -1 1 -1 -1 0];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130301_2
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_2'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 1 -1 1 1 -1 -1 -1 1 1 1 1 1 -1 -1 -1 -1 1 1 -1 1 -1 -1 1 1 0 0 1 -1 1 0 1 -1 -1 0 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130302
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
onOff=[-1 1 0 0 0 0 0 -1 0 0 1 1 -1 1 1 1 -1 1 0 1 -1 1 1 -1 1 1 -1 1 -1 1 -1 1 1 -1 0 -1 1 -1 1];
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130302_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302_1'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 0 0 -1 0 -1 1 -1 1 -1 -1 -1 1 1 -1 1 1 1 -1 0 1 1 -1 0 1 -1 -1 1 -1 0 -1 -1 -1 0 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])



%% 20130220
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 -1 -1 -1 -1 -1 1 1 0 1 -1 -1 1 1 1 1 0 -1 1 1 1 1 1 1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 1 -1 -1 -1 1 0 1 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])



%% 20130220_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 -1 -1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20130224
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130224'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 1 1 1 -1 -1 1 1 1 -1 0 1 -1 -1 1 1 1 0 -1 -1 1 1 -1 1 -1 0 1 1 -1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20130225
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130225'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 1 -1 1 1 1 -1 -1 1 -1 1 1 -1 1 -1 1 0 1 1 -1 -1 1 1 1 1 1 1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

%% 20130226
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130226'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[-1 -1 1 -1 -1 1 -1 -1 -1 0 -1 1 1 -1 1 -1 1 -1 1 -1 -1 -1 1 1 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])



%% 20130227
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130227'
codeWord='HL10'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 1 0 -1 1 1 -1 0 0 0 1 1 -1 -1 0 1 -1 1 0 -1 0 0 0 0 1 1 0 0 -1 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20120329
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120329'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[1 1 -1 1 -1 1 -1 -1 1 -1 0 1 1 0 -1 1 0 0 0 1 -1 1 1 1 -1 0 0 1 0 0 1 1 0 1 -1 -1 1 -1 -1 -1 0 1 1 1 1 1 1 -1 1 1 -1 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20120627
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120627'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[0 -1 1 -1 1 0 1 -1 -1 1 1 1 1 -1 1 1 -1 0];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])


%% 20120714
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120714'
codeWord='ffflicker'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
onOff=[0 1 1 1 0 0 1 1 0 0 0 1 0 1 -1];
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])

