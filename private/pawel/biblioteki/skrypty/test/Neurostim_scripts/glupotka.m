[NumberOfMovies3,NumberOfPatternsPerMovie3,AllPatternsUsed3,Patterns3]=NS512_MoviePatterns('G:\2010-09-11-0\movie003',NS_GlobalConstants);
[NumberOfMovies5,NumberOfPatternsPerMovie5,AllPatternsUsed5,Patterns5]=NS512_MoviePatterns('G:\2010-09-11-0\movie005',NS_GlobalConstants);

Pt=16;

plot([1:512],Patterns3(Pt,:),'bd-',[1:512],Patterns3(Pt,:)-0.1,'gd-')