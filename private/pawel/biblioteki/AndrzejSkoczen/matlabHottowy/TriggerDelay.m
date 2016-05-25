function bitstr = TriggerDelay(Record, Connect, DAC)
            cmd = '100_';
            if(Record < 1023)
                recordstr = dec2bin(Record,10);
            else
                error('error in TriggerDelay: Record trigger time is to long,  should be <= 1023');
            end
            if(Connect < 1023)
                connectstr = dec2bin(Connect,10);
            else
                error('error in TriggerDelay: Connect trigger time is to long,  should be <= 1023');
            end
            if(DAC < 1023)
                DACstr = dec2bin(DAC,10);
            else
                error('error in TriggerDelay: DAC trigger time is to long,  should be <= 1023');
            end
            bitstr = ['1010_',cmd,recordstr,'_',connectstr,'_',DACstr];
end
