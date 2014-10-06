%cd D:\analysis\2008-08-27-4\data002;
cd /Volumes/Lee/Analysis/Lauren/2008-08-26-0;
fid=fopen('clusters006artifacts','wb','ieee-le');
NumberOfMovies=76;
NumberOfPatterns=64;
NumberOfRepetitions=100;


a=ones(NumberOfMovies*NumberOfPatterns*NumberOfRepetitions+3,1);
a(1)=NumberOfMovies; %126
a(2)=NumberOfPatterns; %330
a(3)=NumberOfRepetitions;
fwrite(fid,a,'int16');
fclose(fid)