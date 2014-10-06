clear;

kanal1=struct('kanal',1,'data',rand(5,50));
kanal2=struct('kanal',2,'data',rand(5,60));
kanal3=struct('kanal',3,'data',rand(5,70));

pattern1=[kanal1 kanal2 kanal3];
pattern2=[kanal2 kanal2];

chunk=[pattern1 pattern2]