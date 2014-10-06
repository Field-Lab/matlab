function f=szumy_filtr(f1,f2,f);

f0=(1+f/f1).^(-4);
f1=f0.*(f/f2).^4;
f2=f1.*(1+f/f2).^(-4);
f=f2./max(max(f2));
