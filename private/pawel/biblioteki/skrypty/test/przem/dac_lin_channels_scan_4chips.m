clear

cd H:/pliki/nauka/stymulacja/chip/testy/2006_07_12/

a=importdata('60nA_4700k_channels_scan_board1_chip_left.dat');
b=importdata('4uA_75k_channels_scan_board1_chip_left.dat');
b=importdata('250uA_1190_channels_scan_board1_chip_left.dat');
slope1=przekladki_channels_scan(a,b);

a=importdata('60nA_4700k_channels_scan_board1_chip_right.dat');
b=importdata('4uA_75k_channels_scan_board1_chip_right.dat');
slope2=przekladki_channels_scan(a,b);

a=importdata('60nA_4700k_channels_scan_board3_chip_left.dat');
b=importdata('4uA_75k_channels_scan_board3_chip_left.dat');
b=importdata('250uA_1190_channels_scan_board3_chip_left.dat');
slope3=przekladki_channels_scan(a,b);

a=importdata('60nA_4700k_channels_scan_board3_chip_right.dat');
b=importdata('4uA_75k_channels_scan_board3_chip_right.dat');
slope4=przekladki_channels_scan(a,b);

figure(11)
s=[1:32];
plot(s,slope1,'bd-',s,slope2,'kd-',s,slope3,'rd-',s,slope4,'cd-');