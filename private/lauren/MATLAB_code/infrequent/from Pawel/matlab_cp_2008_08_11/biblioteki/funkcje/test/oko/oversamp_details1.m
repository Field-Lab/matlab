function y=oversamp_details1(s,figura,marg_left,marg_right);
%swgtwrgt
prog=200;
histereza=40;
znak=-1;

N=2;
filtr_dl=16;
f_gr=0.98;

dl=length(s);
signal=s(1:2:dl); 

%3. Oversampling:
[y,filtr]=oversampling(signal,N,filtr_dl,f_gr);
y1=y(1,filtr_dl+1:filtr_dl+dl);
roznica=y1-s;
    
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
wynik0=detect2(s(marg_left:dl-marg_right),detect_param); 
wynik0=wynik0+marg_left; %bo oczywiste

l_wynik0=length(wynik0);
for i=1:l_wynik0
    %subplot(1,2,1);
    %plot(s(wynik0(1,i)-marg_left:wynik0(1,i)+marg_right));
    %hold on;
    %subplot(1,2,2);
    loglog(abs(fft(s(wynik0(1,i)-marg_left:wynik0(1,i)+marg_right))));
    axis([0 round((marg_left+marg_right)/2) 1 20000]);
    hold on;
end
y=y1;