function bitstr=ChannelMask_PH(nrch,DL,CM);

DLstr = dec2bin(DL,12);
bitstr=['1010_' '011_' DLstr '_' RejectSpaces(num2str(CM))];