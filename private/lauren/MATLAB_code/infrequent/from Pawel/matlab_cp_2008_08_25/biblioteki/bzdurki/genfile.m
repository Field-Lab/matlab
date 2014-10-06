cd /Volumes/Lee/Analysis/Lauren/2008-08-27-4;
fid=fopen('clusters007','wb','ieee-le');
NumberOfMovies=26;
NumberOfPatterns=124;
NumberOfRepetitions=100;
a=ones(NumberOfMovies*NumberOfPatterns*NumberOfRepetitions+3,1);
a(1)=NumberOfMovies; %126
a(2)=NumberOfPatterns; %330
a(3)=NumberOfRepetitions;
fwrite(fid,a,'int16');
fclose(fid)