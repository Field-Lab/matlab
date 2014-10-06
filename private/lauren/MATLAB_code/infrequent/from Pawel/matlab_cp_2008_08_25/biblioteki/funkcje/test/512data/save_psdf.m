function w=save_psdf(source_path,destin_path,filename,channel0,f0,fs,wykladnik,length);

source_path
cd(source_path);
%outname=[filename(1:7) 'psdf' '.bin']
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik,length);
data=[abs(y') sigma' psdf];

cd(destin_path);
outname=[filename(1:7) 'psdf' '.bin']
w=outname;

fid=fopen(outname,'wb')
fwrite(fid,data,'double');
fclose(fid);