function Gain=NS256_LoadGainValue(Filename);

GainValues=[160 270 360 440 510 560 610 650 680 720 750 760 770 790 800 810 820 830 835 840];

c1=fopen(Filename','r');
a1=fread(c1);
fclose(c1);

b1=abs((a1(145:149)-1)');
GainCode=b1(1)*16+b1(2)*8+b1(3)*4+b1(4)*2+b1(5)
if GainCode>length(GainValues)
    error('GainCode too large');
    Gain=0;
else
    Gain=GainValues(GainCode);
end