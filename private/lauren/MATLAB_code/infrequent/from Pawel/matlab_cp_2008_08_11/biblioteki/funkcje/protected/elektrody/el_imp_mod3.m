function imp=el_imp_mod3(Ce,beta,Re,Rs,Cpar,f);

zcpe=Ce*(i*2*pi*f).^(-beta);
a=zcpe+Re;
b=zcpe*Re;
z1=b./a;
imp0=z1+Rs;

Zpar=(Cpar*i*2*pi*f).^(-1);
a=imp0+Zpar;
b=imp0.*Zpar;
imp=b./a;
%imp=Zpar;