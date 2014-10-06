
cd H:/pliki/nauka/stymulacja/chip/testy/CURRENT_MODE/

dane0=importdata('curr_r1_19200k_b4_c1.dat')';
dac=[-127:127];
R1=19200e3;

doubleaxis(dac,dane0,R1);
