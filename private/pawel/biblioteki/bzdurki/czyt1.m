path1='G:\analysis\slices\2011-06-29-1\data001\data001000.bin'
fid = fopen(path1, 'r+');
offset = 2200000*770
    %if FileNumber==0
    %    offset=offset+size(binRep,1);
    %end
fseek(fid, offset, 'bof');
a=fread(fid, 770*10000,'uint8');
fclose(fid);