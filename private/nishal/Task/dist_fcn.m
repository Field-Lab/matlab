function z = dist_fcn(ji,j_cuei,t,J_card,Z_card)
 
zz = dec2base(t,Z_card,J_card);
zz=zz(end-J_card+1:end);
zz = str2double(zz(j_cuei))-1;

ztrue = sign(ji-j_cuei);

z=double(zz~=ztrue);
end