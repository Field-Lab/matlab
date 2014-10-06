nazwa_chipu='chipB_';

g_name=[nazwa_chipu '_gain.dat']
o_name=[nazwa_chipu '_offset.dat']
k_name=[nazwa_chipu '_korr.dat']
s_name=[nazwa_chipu '_setup.dat']

gain=importdata(g_name);
offset=importdata(o_name);
korr=importdata(k_name);
setup=importdata(s_name);

tlumik=setup(1,1);
sl=length(setup);
gen_ampl=setup(1,2:sl);

sgain=size(gain);
for i=1:sgain(2)
  gain(:,i)=gain(:,i)./gen_ampl(i);
end
gain=gain.*tlumik;
  
a=max(gain');
kanaly=find(a>100);

wzm_fit=zeros(length(kanaly),2)

pocz=5; %od ktorej amplitudy wej. bierze sie do fitu
koniec=14; % do ktorej

for i=1:length(kanaly)
  wzm=gain(kanaly(i),pocz:koniec)
  [pg,sg]=polyfit(gen_ampl(1,pocz:koniec),wzm,1);
  wzm_fit(i,1:2)=pg;
  a0(i)=sum(wzm)/sum(gen_ampl(1,pocz:koniec));
end

nachylenie=wzm_fit(:,2);
dobre=find(nachylenie>100);
figure(21)
plot(wzm_fit(:,2));

figure(52)
plot(a0)

mean(wzm_fit(dobre,2))
std(wzm_fit(dobre,2))

mean(a0)
std(a0)





