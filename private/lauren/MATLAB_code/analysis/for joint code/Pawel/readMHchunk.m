function header=readMHchunk(fid0)
%Reads the header-of-movie-file chunk and returns the two-field structure:
%header=struct(number_of_movies,period), where 'number_of_movies' is the
%total number of Movie Data (MD) chunks in the file and the 'period' value
%is the period of repetition for each movie in the file. If the 'period'
%value equals to 0, then the period is defined individually for each movie
%by the value saved in the MD chunk.

%fid0=fopen(filename,'r','b');
ID=fread(fid0,8,'int8')';

if (ID~=[86 17 -34 51 11 -15 49 30])
   error('chunk ID not correct');
end

header=ID;
size=fread(fid0,1,'int64');
number_of_movies=fread(fid0,1,'int32');
period=fread(fid0,1,'int32');
header=struct('number_of_movies',number_of_movies,'period',period);

