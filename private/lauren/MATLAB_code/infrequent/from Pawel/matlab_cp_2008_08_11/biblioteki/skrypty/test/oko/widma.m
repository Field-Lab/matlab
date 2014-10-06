fp=20000;
name='Data008conv';
header=206;
nrchns=65;

channels1=[2:17];
channels2=channels1+16;
channels3=channels1+32;
channels4=channels1+48;
samples=[1 20000];

r1=readconv(name,header,nrchns,channels1,samples);
r2=readconv(name,header,nrchns,channels2,samples);
r3=readconv(name,header,nrchns,channels3,samples);
r4=readconv(name,header,nrchns,channels4,samples);

lsb=0.001;
N=4096;
gwsz=zeros(16,N);

for i=1:16
   s=r2(i,:);
   srednia=mean(s)
   [f,w]=spdf(s,N,N,fp,lsb);
   gwsz(i,:)=w;
end
