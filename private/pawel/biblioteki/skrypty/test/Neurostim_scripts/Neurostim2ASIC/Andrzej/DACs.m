classdef DACs
    enumeration
        Ipol1, Ipol2, Flow, Fhigh, Gain, n_u_1, n_u_2, Istim
    end
    methods
        function bitstr = AnalogSettings(PA, PV)
            cmd = '000_';
            switch(PA)
                case(DACs.Ipol1)
                    pastr = '000';
                case(DACs.Ipol2)
                    pastr = '001';
                case(DACs.Flow)
                    pastr = '010';
                case(DACs.Fhigh)
                    pastr = '011';
                case(DACs.Gain)
                    pastr = '100';
                case(DACs.Istim)
                    pastr = '111';
            end
            if(PV < 1023)
                pvstr = dec2bin(PV,10);
            else
                error('error in AnalogSettings: PV is too big, should be <= 1023');
            end
            bitstr = ['1010_',cmd,pastr,'_',pvstr];
        end
    end
end
