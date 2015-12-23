function d = dist_fcn_2(ji,j_cuei,t,Z_card,Dict,P,Q)
 
I=element_in_Dict(ji,j_cuei,Dict);

ztrue=I;

z_dec=t_fcn(t,j_cuei,Z_card,Q);

d = double(ztrue~=z_dec);
end

function z_dec = t_fcn(t,j_cuei,Z_card,Q) 
%No of fcns = Z_Card^Q(=J_card??)
zz = dec2base(t,Z_card,Q);
zz=zz(end-Q+1:end);
z_dec = str2double(zz(j_cuei));

end