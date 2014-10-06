 cd /home/pawel/pliki/nauka/neuro64/
 cd lipiec2002

 names(1)={'dane_20khz.dat'};
 names(2)={'dane_22khz.dat'};
 names(3)={'dane_23khz.dat'};
 names(4)={'dane_25khz.dat'};
 names(5)={'dane_28khz.dat'};
 names(6)={'dane_30khz.dat'};
 %n=[40 40 400 40];
 n=40;
 fp=[20000 22000 23000 25000 28000 30000];
 
 N=10000;
 %dane=zeros(4,16384,64);
 
 for i=1:6
     name=names{i}
     a=importdata(name);
     s=a(1:round(N/2),:)'*4;  %mnozenie przez 4 - kompensacja Labviewwowj f-cji do gest. widm.
     s=s/n;  %40 - ilosc powtorzen pomiarow
 %s=s/N*fp(i);
 
     f=[1:round(N/2)]/N*fp(i);

     signal=s(35,:);
     subplot(3,2,i);
     %loglog(f,signal);
     semilogy(f,signal);
     axis([min(f) max(f) 1e-10 1e-5]);
     h=gca;
     set(h,'FontSize',12);
     grid on;
%plot(f,signal);
     %sqrt(sum(signal))
%axis([1 10000 1e-12 1e-2]);
end
