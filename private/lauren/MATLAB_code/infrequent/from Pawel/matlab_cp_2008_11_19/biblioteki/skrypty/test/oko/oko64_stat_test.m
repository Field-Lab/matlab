cd /home/pawel/pliki/oko/oko_64/koncowe/17-12-2001;
tlumik=2.12/0.026
gen_ampl=0.03;

figure(1);
%clf;
figure(2);
%clf;

nrchns=64;
name_dane='chipB_gr4_data.dat';
name_clk='chipB_gr4_clk.dat';

%name_dane='chipA_gr3_dane_highres.dat';
%name_clk='chipA_gr3_clk_highres.dat';

a=importdata(name_dane);
c=importdata(name_clk);

l=length(a);
signal=a(5:l,1)';
clock=c(5:l,1)';

w=oko64_stat(signal,clock,nrchns);
prog=0.3;
A=0.4;
f=0.039;
fi=1;
offset=-0.2;
b=[A f fi offset];

ch_tekst='channel: ';
gain_tekst='gain: ';
offset_tekst='offset: ';
czest_tekst='czest: ';
ch=0;
figure(1);
for i=1:nrchns
    %subplot(8,8,i);
    s=w(i,:);
    t=[1:length(s)];
    ampl=max(s)-min(s);
    if ampl>prog
        i;
        [b,p]=nielfit(t,s,'sinus_fit',b,[0.0005 0.0002 0.0005 0.0002]);
        b(1,4);
        aproks=feval('sinus_fit',b,t);
        
        figure(1);
        subplot(8,8,i);
        plot(t,s,t,aproks);
        axis([1 32 -.5 .5]);
        grid on;
        
        figure(2);
        ch=ch+1;
        subplot(5,5,il);
        il=il+1;
        plot(t,s,t,aproks);
        axis([1 32 -.5 .5]);
        
        tekst=num2str(i);
        text(4,0.4,[ch_tekst tekst]);
        tekst=num2str(b(1,1)/gen_ampl*tlumik,'%2.3f');
        text(4,0.25,[gain_tekst tekst]);
        tekst=num2str(b(1,4),'%2.3f');
        text(4,0.1,[offset_tekst tekst ' mV']);
        
        tekst=num2str(b(1,2),'%2.4f');
        %text(4,0.1,[czest_tekst tekst]);
        
        grid on;
    else 
        figure(1);
        subplot(8,8,i);
        plot(w(i,:));
        axis([1 32 -.5 .5]);
        grid on;
    end
    %axis([1 32 -.5 .5]);
end


