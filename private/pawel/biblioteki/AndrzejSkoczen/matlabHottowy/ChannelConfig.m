function bitstr = ChannelConfig(nrch,CHCL)
%nrch - liczba kanalow w modelu, ktory zamierzamy symulowac
%dla pelnego modelu nrch=64 i mozna stowac wszystkie nizsze wartosci
    CHC = '';
    maxCHCL = nrch*14;
    cmd = '010_';
    if(CHCL <= maxCHCL)
        CHCLstr = dec2bin(CHCL,10);
    else
        error('error in ChannelConfig: CHCL is too big,  should be <= %i',maxCHCL);
    end
    for i = 1:CHCL
       CHC = [CHC,'1']; 
    end
    bitstr = ['1010_',cmd,CHCLstr,'_',CHC];
end