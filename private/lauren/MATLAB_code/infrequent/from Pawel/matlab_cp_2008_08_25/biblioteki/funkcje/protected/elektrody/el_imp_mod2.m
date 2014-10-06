function imp=el_imp_mod2(Ce1,beta1,Re,Ce2,beta2,f);

zcpe1=Ce1*(i*2*pi*f).^(-beta1);
zcpe2=Ce2*(i*2*pi*f).^(-beta2);
a=zcpe1+Re;
b=zcpe1*Re;
z1=b./a;
imp=z1+zcpe2;