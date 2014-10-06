function [V,Icpe,Ire]=art_cpe_nielin_sym1_triphas(model_nielin,A1,t1,A2,t2,A3,t3,T,dt);

[V1,Icpe1,Ire1]=art_cpe_nielin_sym1(model_nielin,A1,t1,T,dt,0);
'ble'
Vpocz=V1(floor(t1/dt));
[V2,Icpe2,Ire2]=art_cpe_nielin_sym1(model_nielin,A2,t2,T-t1,dt,Vpocz);
'bleble'
Vpocz=V2(floor(t2/dt));
[V3,Icpe3,Ire3]=art_cpe_nielin_sym1(model_nielin,A3,t3,T-t2-t1,dt,Vpocz);

V2=[zeros(1,length(V1)-length(V2)) V2];
V3=[zeros(1,length(V1)-length(V3)) V3];
V=V1+V2+V3;

Icpe2=[zeros(1,length(Icpe1)-length(Icpe2)) Icpe2];
Icpe3=[zeros(1,length(Icpe1)-length(Icpe3)) Icpe3];
Icpe=Icpe1+Icpe2+Icpe3;

Ire2=[zeros(1,length(Ire1)-length(Ire2)) Ire2];
Ire3=[zeros(1,length(Ire1)-length(Ire3)) Ire3];
Ire=Ire1+Ire2+Ire3;