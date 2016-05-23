function frame=insert_current(channel, frame, amp, range, r_amp, amp_max,number)

%connect
%hold switch3 zawsze na 1
%discharge switch2 zawsze na 0
%record switch1
%state=[connect;switch1;switch3;switch2];
%number - numer ramki
% state =[0;1;1;0];
%         for i=1:64
%             frame(:,(3*i-2)+(250*(number-1)))=state;
%         end
        
prad=dac(amp, range, r_amp, amp_max);
p=prad(:,1);
d=prad(:,2);
frame(:,((3*channel-1)+250*(number)):(3*channel+250*(number)))=[p,d];






    
        
    



  


