function sample=DAC_and_switches_decode(DAC_and_switches)
%ble ble

switches=floor(DAC_and_switches/65536);
DAC=DAC_and_switches-(switches*65536)-128;
discharge=switches-2*floor(switches/2);
switches=(switches-discharge)/2;
hold=switches-2*floor(switches/2);
switches=(switches-hold)/2;
record=switches-2*floor(switches/2);
switches=(switches-record)/2;
connect=switches-2*floor(switches/2);
%DAC=DAC_and_switches-(switches*65536)-128

sample=struct('DAC',DAC,'record',record,'connect',connect,'discharge',discharge,'hold',hold);