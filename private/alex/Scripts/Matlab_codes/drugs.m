%saving
for i=1
OnSpikes08=OffSpikes;
OffSpikes08=OffSpikes;
Spont08=spontSpikes;
times08=[1,19,54,67,85,109,159];


OnSpikes11=OnSpikes;
OffSpikes11=OffSpikes;
Spont11=spontSpikes;
times11=[1,14,46,90,106,178,234];

OnSpikes15=OnSpikes;
OffSpikes15=OffSpikes;
Spont15=spontSpikes;
times15=[1,23,71,121,177,263,320];

OnSpikes19=OnSpikes;
OffSpikes19=OffSpikes;
Spont19=spontSpikes;
times19=[1 26 68 106 151 276 282];

OnSpikes16=OnSpikes;
OffSpikes16=OffSpikes;
Spont16=spontSpikes;
times16=[1 42 79 128 176 225 271 319 360 399 420 440 482 517 533 563 614 649];

save('C:\Documents and Settings\atikidzhi\Desktop\all_on_off','OnSpikes08','OffSpikes08','Spont08','OnSpikes11','OffSpikes11','Spont11','OnSpikes15','OffSpikes15','Spont15','OnSpikes16','OffSpikes16','Spont16','OnSpikes19','OffSpikes19','Spont19','times08','times11','times15','times16','times19')
end

load('C:\Documents and Settings\atikidzhi\Desktop\all_on_off.mat')

a=cell(1,15);
a{1}='OnSpikes08';
a{2}='OnSpikes11';
a{3}='OnSpikes15';
a{4}='OnSpikes16';
a{5}='OnSpikes19';
a{6}='OffSpikes08';
a{7}='OffSpikes11';
a{8}='OffSpikes15';
a{9}='OffSpikes16';
a{10}='OffSpikes19';
a{11}='Spont08';
a{12}='Spont11';
a{13}='Spont15';
a{14}='Spont16';
a{15}='Spont19';


contr_08=mean(OnSpikes08(2:15)-Spont08(2:15))
APB10_08=mean(OnSpikes08(25:50)-Spont08(25:50))
APB10_2nd_08=mean(OnSpikes08(69:83)-Spont08(69:83))
wash_08=mean(OnSpikes08(90:105)-Spont08(90:105))

contr_11=mean(OnSpikes11(2:12)-Spont11(2:12))
APB10_11=mean(OnSpikes11(25:44)-Spont11(25:44))
APB10_2nd_11=mean(OnSpikes11(93:103)-Spont11(93:103))
wash_11=mean(OnSpikes11(130:170)-Spont11(130:170))

contr_15=mean(OnSpikes15(2:22)-Spont15(2:22))
CPP10_15=mean(OnSpikes15(44:69)-Spont15(44:69))
CPP10_2nd_15=mean(OnSpikes15(130:170)-Spont15(130:170))
wash_15=mean(OnSpikes15(225:260)-Spont15(225:260))

contr_16=mean(OnSpikes16([3:40,100:125])-Spont16([3:40,100:125]))
CPP10_16=mean(OnSpikes16(155:174)-Spont16(155:174))
CPP10_2nd_16=mean(OnSpikes16(250:269)-Spont16(250:269))
wash1_16=mean(OnSpikes16([300:318,380:400])-Spont16([300:318,380:400]))
CPP20_16=mean(OnSpikes16(410:419)-Spont16(410:419))
CPP20_2nd_16=mean(OnSpikes16(455:475)-Spont16(455:475))
ALL20_16=mean(OnSpikes16(495:515)-Spont16(495:515))
ALL20_2nd_16=mean(OnSpikes16(538:562)-Spont16(538:562))
wash_16=mean(OnSpikes16([600:612,625:649])-Spont16([600:612,625:649]))

contr_19=mean(OnSpikes19(50:65)-Spont19(50:65))
ALL20_19=mean(OnSpikes19(80:100)-Spont19(80:100))
ALL100_19=mean(OnSpikes19(120:140)-Spont19(120:140))
wash_19=mean(OnSpikes19(220:250)-Spont19(220:250))


contr_08=mean(OffSpikes08(2:15)-Spont08(2:15))
APB10_08=mean(OffSpikes08(25:50)-Spont08(25:50))
APB10_2nd_08=mean(OffSpikes08(69:83)-Spont08(69:83))
wash_08=mean(OffSpikes08(90:105)-Spont08(90:105))

contr_11=mean(OffSpikes11(2:12)-Spont11(2:12))
APB10_11=mean(OffSpikes11(25:44)-Spont11(25:44))
APB10_2nd_11=mean(OffSpikes11(93:103)-Spont11(93:103))
wash_11=mean(OffSpikes11(130:170)-Spont11(130:170))

contr_15=mean(OffSpikes15(2:22)-Spont15(2:22))
CPP10_15=mean(OffSpikes15(44:69)-Spont15(44:69))
CPP10_2nd_15=mean(OffSpikes15(130:170)-Spont15(130:170))
wash_15=mean(OffSpikes15(225:260)-Spont15(225:260))

contr_16=mean(OffSpikes16([3:40,100:125])-Spont16([3:40,100:125]))
CPP10_16=mean(OffSpikes16(155:174)-Spont16(155:174))
CPP10_2nd_16=mean(OffSpikes16(250:269)-Spont16(250:269))
wash1_16=mean(OffSpikes16([300:318,380:400])-Spont16([300:318,380:400]))
CPP20_16=mean(OffSpikes16(410:419)-Spont16(410:419))
CPP20_2nd_16=mean(OffSpikes16(455:475)-Spont16(455:475))
ALL20_16=mean(OffSpikes16(495:515)-Spont16(495:515))
ALL20_2nd_16=mean(OffSpikes16(538:562)-Spont16(538:562))
wash_16=mean(OffSpikes16([600:612,625:649])-Spont16([600:612,625:649]))

contr_19=mean(OffSpikes19(50:65)-Spont19(50:65))
ALL20_19=mean(OffSpikes19(80:100)-Spont19(80:100))
ALL100_19=mean(OffSpikes19(120:140)-Spont19(120:140))
wash_19=mean(OffSpikes19(220:250)-Spont19(220:250))
