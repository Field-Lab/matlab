function y=xls2bin(filename_r,filename_w,blocksize);

i=0;
fid=fopen(filename_w,'w');
s=blocksize;
while s==blocksize
	i=i+1;
	a=readxls(filename_r,(i-1)*blocksize,blocksize);
	s=length(a);
	fwrite(fid,a,'double');
end
fclose(fid);
