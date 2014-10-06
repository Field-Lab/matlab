NS_GlobalConstants = NS_GenerateGlobalConstants(61);
MovieNumber=42;

DataMovieFile='C:\pawel\nauka\Stimchip_paper\2012grudzien\dyskusje\2012_01_05\movie017'
MovieData1=NS512_MovieData2(DataMovieFile,MovieNumber,NS_GlobalConstants);

a=find(MovieData1==61);

for i=1:length(a)
    time=MovieData1(a(i)-1);
    b=find(MovieData1==time | MovieData1==time-40 | MovieData1==time-80 | MovieData1==time-120 | MovieData1==time-160);
    if length(b)==1
        i
    end
end

%only el 61: pulse 2,3,5,7,8,13,16,20,21; excluding cases when other
%electrode stimulated 2 ms earlier: 3,5,7,8,21; excluding difference of
%40ms: 3,7,8,21; 7,8,21 - wiecej niz 8 ms;