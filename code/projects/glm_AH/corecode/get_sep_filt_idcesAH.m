function [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idcesAH(offset,Basepars)

[n mk] = get_nspacetimeAH(Basepars);
        


s1_idx = offset+1:offset+n;
offset = offset + n;
t1_idx = offset+1:offset+mk;
offset = offset + mk;
s2_idx = offset+1:offset+n;
offset = offset + n;
t2_idx = offset+1:offset+mk;

