function ThresholdFile=NS512_GenerateExampleThresholdFile(Path);

a=zeros(12,4);
a(:,1)=[76 227 256 271 391 406 541 616 691 736 856 901];
a(:,2)=[129:10:239];
a(:,3)=a(:,2);
a(:,4)=rand(1,12)+0.2;

fid1=fopen(Path,'w');
fwrite(fid1,a,'double');
ThresholdFile=fclose(fid1);