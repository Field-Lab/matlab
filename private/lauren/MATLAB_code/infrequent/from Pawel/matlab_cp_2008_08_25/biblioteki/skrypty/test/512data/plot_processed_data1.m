cd H:/pliki/nauka/crosstalk/;

for i=1:40
    if (i<10)
        przedrostek='data00';
    else
        przedrostek='data0';
    end
    filename=[przedrostek num2str(i) '.bin'];
    fid=fopen(filename,'r');
    fseek(fid,512*16384*8,-1);
    h=fread(fid,[512 2],'double');
    fclose(fid);
    
    figure(11);
    subplot(4,10,i);
    chns=[1:512];
    semilogy(chns,abs(h(:,1)),'b-',chns,h(:,2),'r-');
%    semilogy(chns,abs(h(:,1))')
    axis([0 512 1e-2 1e3]);
    grid on;            
end

filename='data029.bin';
f=[0:16383]/16384*20000;
figure(12);
fid=fopen(filename,'r');
%fseek(fid,300*16386*8,-1);
h=fread(fid,[512 16386],'double');       
   for i=1:25
       subplot(5,5,i);
       loglog(f,h(i,1:16384));
       axis([1 10000 1e-4 10]);
       grid on;
   end
figure(13);
fseek(fid,512*16384*8,-1);
h=fread(fid,[512 2],'double');
plot(chns,abs(h(:,1)),'bd-',chns,h(:,2),'rd-');
axis([0 512 0.01 1000]);
grid on;
fclose(fid);    

filename='data051.bin';
f=[0:16383]/16384*20000;
figure(14);
fid=fopen(filename,'r');
%fseek(fid,300*16386*8,-1);
h=fread(fid,[512 16386],'double');       
   for i=1:25
       subplot(5,5,i);
       loglog(f,h(i+140,1:16384));
       axis([1 10000 1e-4 100]);
       grid on;
   end
   
figure(15);
fseek(fid,512*16384*8,-1);
h=fread(fid,[512 2],'double');
plot(chns,h(:,2),'rd-');
axis([0 512 0.01 70]);
grid on;
fclose(fid);    