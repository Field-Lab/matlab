fid = fopen('C:\home\Pawel\nauka\CMOS\Neurostim3\DataStreamAnalogSettings.txt','w');
HW = 40;
HD = 0;
nrch = 64;

ParameterID=0;
bitstr = []; %[SoftReset, '__',HoldTiming(HW,HD),'__',TriggerDelay(10,5,6),'__'];
for i=1:10
    bitstr=[bitstr AnalogSettingsPH(ParameterID,2^(i-1))];
end

for i=1:7
    bitstr=[bitstr AnalogSettingsPH(i,255)];
end

bitstr   
fprintf(fid,'%s',bitstr);
fclose(fid);