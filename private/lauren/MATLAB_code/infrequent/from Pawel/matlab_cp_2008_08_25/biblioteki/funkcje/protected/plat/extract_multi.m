function dane=extract_multi(p,numer,nrchns);

size(p);

start=(numer-1)*4; %bo 4 kolumny na jeden pomiar

data0a=p(:,start+2);
data0b=p(:,start+3);
data1=data0a(5:length(data0a),1)'-data0b(5:length(data0a),1)';
clear data0a;
clear data0b;

clk0=p(:,start+4);
clk1=clk0(5:length(clk0),1)';
clear clk0;

trig0=p(:,start+1);
trig=trig0(5:length(trig0),1)';
clear trig0;

%trig=trig-(max(trig)-min(trig))/2;
trig=trig-2;   
s_trig=sign(trig);
sd_trig=diff(s_trig); %badamy, w ktorym miejscu sygnal trigera przekroczyl prog


'hgvkv';
size(sd_trig);
a=find(sd_trig>0);
size(a);
if length(a)~=2
  warning('klopoty z trigerem');
  a0=a;
  clear a;
  a(1,1)=a0(1,1);
  a(1,2)=a0(1,length(a0));
end
size(a);

'ggggg';


clk=clk1(1,a(1,1):a(1,2));
clear clk1;
size(clk);
clk=clk-1.2;
%clk=clk-(max(clk)-min(clk))/2;
%figure(3)
%plot(clk)

data=data1(1,a(1,1):a(1,2));
%plot(data)
clear data1;
size(data);

s_clk=sign(clk);
s=s_clk+abs(s_clk);
ds=diff(s);
size(ds);
p=find(ds==2);
size(p);

dlugosc=floor(length(p)/nrchns);

dane=zeros(nrchns,dlugosc);
for i=1:nrchns
  chns(i,1:dlugosc)=data(p(i:nrchns:dlugosc*nrchns));
end

dane=chns;

