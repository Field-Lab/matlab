gen_ampl=[0.01:0.004:0.08];

cd /home/pawel/pliki/oko/oko_64/koncowe
name(1)={'czest_chipA_gr2.dat'};
%name(2)={'chipA_gr2_370Hz_2.dat'};
%name(3)={'chipA_gr3_370Hz_2.dat'};
%name(4)={'chipA_gr4_370Hz_2.dat'};
sname=size(name);

%d=zeros(

for plik=1:sname(2)
  nazwa=name{plik}
  dane=importdata(nazwa);
  %for i=1:length(gen_ampl)
    d=extract_multi(dane,5,64);
  %end
end

figure(1)
for i=1:64
  subplot(8,8,i)
  plot(d(i,:))
end
