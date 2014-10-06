ch1=struct('channel',1,'data',rand(5,10));
ch2=struct('channel',2,'data',rand(5,12));

pattern1=[ch1 ch2];

ch3=struct('channel',3,'data',rand(5,8));
ch4=struct('channel',4,'data',rand(5,12));
ch5=struct('channel',5,'data',rand(5,15));

pattern2=[ch3 ch4 ch5];

patterns=[pattern1 pattern2]