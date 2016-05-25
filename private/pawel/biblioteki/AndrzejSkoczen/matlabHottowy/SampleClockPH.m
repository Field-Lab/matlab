function bitstr = SampleClockPH(CM_DIV,CM_DEL,HA_DEL,HA_INV,HA_W);
%nrch - liczba kanalow w modelu, ktorzy zamierzamy symulowac
%dla pelnego modelu nrch=64 i mozna stowac wszystki inizesz wartosci
    
cmd = '101_';
bitstr=[cmd dec2bin(CM_DIV,10) '_' dec2bin(CM_DEL,10) '_' dec2bin(HA_DEL,10) '_' dec2bin(HA_INV,7) '_' dec2bin(HA_W,10) '_'];