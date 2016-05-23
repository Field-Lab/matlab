function bitstr = HoldTiming(HW, HD)
            cmd = '110_';
            HWstr = dec2bin(HW,6);
            HDstr = dec2bin(HD,10);
            bitstr = ['1010_',cmd,HWstr,'_',HDstr];
end
