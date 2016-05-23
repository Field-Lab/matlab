function [frame_zero]=frame_zero()
hold=1;
discharge=0;
record=1;
connect=0;
state =[connect;record;hold;discharge];
frame_zero=[repmat([state, zeros(4,2)],1,64),zeros(4,58)];
end