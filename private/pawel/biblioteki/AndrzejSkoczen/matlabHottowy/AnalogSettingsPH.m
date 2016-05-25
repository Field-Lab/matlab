function bitstr=AnalogSettingsPH(ParameterID,Value);

cmd='000_';
ParameterAddress=dec2bin(ParameterID,3);
ParameterValue=dec2bin(Value,10);

bitstr=[cmd ParameterAddress '_' ParameterValue '_'];