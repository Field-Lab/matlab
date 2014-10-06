name1='chipB_gr1_all_czest.dat';
name2='chipB_gr2_all_czest.dat';
name3='chipB_gr3_all_czest.dat';
name4='chipB_gr4_all_czest.dat';

%a1=importdata(name1);
%a2=importdata(name2);
%a3=importdata(name3);
%a4=importdata(name4);

for i=1:20
  figure(101);
  subplot(4,5,i);
  plot(a1(:,(i-1)*4+2));
  
  figure(102);
  subplot(4,5,i);
  plot(a2(:,(i-1)*4+2));
  
  figure(103);
  subplot(4,5,i);
  plot(a3(:,(i-1)*4+2));
  
  figure(104);
  subplot(4,5,i);
  plot(a4(:,(i-1)*4+2));
  
end
 
   
