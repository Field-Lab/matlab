cd /home2/pawel/usa_2004/data/March_19_Monkey;

blocksize=100000;
tic
filename_r='read20_stim19_at_threshold.xls';
filename_w='read20_stim19_at_threshold.dat';
y=xls2bin(filename_r,filename_w,blocksize);
toc

filename_r='no_cadmium_rec12_stim15.xls';
filename_w='no_cadmium_rec12_stim15.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='no_cadmium_rec13_stim17.xls';
filename_w='no_cadmium_rec13_stim17.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='no_cadmium_rec15_stim41.xls';
filename_w='no_cadmium_rec15_stim41.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='no_cadmium_rec20_stim19.xls';
filename_w='no_cadmium_rec20_stim19.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='no_cadmium_rec64_stim49.xls';
filename_w='no_cadmium_rec64_stim49.dat';
y=xls2bin(filename_r,filename_w,blocksize);


filename_r='with_cadmium_rec12_stim15.xls';
filename_w='with_cadmium_rec12_stim15.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='with_cadmium_rec13_stim17.xls';
filename_w='with_cadmium_rec13_stim17.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='with_cadmium_rec15_stim41.xls';
filename_w='with_cadmium_rec15_stim41.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='with_cadmium_rec20_stim19.xls';
filename_w='with_cadmium_rec20_stim19.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='with_cadmium_rec64_stim49.xls';
filename_w='with_cadmium_rec64_stim49.dat';
y=xls2bin(filename_r,filename_w,blocksize);

cd 5_piece;

filename_r='read15_stim22.xls';
filename_w='read15_stim22.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='read15_stim41.xls';
filename_w='read15_stim41.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='read31_stim32_2.xls';
filename_w='read31_stim32_2.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='read31.xls';
filename_w='read31.dat';
y=xls2bin(filename_r,filename_w,blocksize);

filename_r='read_rozne.xls';
filename_w='read_rozne.dat';
y=xls2bin(filename_r,filename_w,blocksize);

toc
