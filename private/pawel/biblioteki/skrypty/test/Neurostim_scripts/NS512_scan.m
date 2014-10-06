[electrodes,Array,MovieChunksFile]=NS512_1el_scan_simple(10,7.5,10000);

cd C:\home\pawel\praca\stim_files; 
fid = fopen('NS512_1el_simple_el','wb');
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('NS512_1el_simple_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('NS512_1el_simple_mv','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 