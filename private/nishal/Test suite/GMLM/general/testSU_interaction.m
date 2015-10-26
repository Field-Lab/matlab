
function testSU_interaction(u1,u2,mask,mov,Y_C,jcell,probe,thperc,pts,type,colr)


up1 = zeros(sum(mask),1);up1(u1)=1;
 up2 = zeros(sum(mask),1);up2(u2)=1;
 inp1 = up1'*mov(:,probe);
 inp2 = up2'*mov(:,probe);
 resp = Y_C(jcell,probe);
 thip2 = prctile(inp2,thperc);
%  pts=7;
if(strcmp(type,'geq'))
 idxx = inp2>=thip2;
else
 idxx=inp2<=thip2;
end

 aa=inp1(idxx);
  plotio_curve_meancorrected(aa,resp(idxx),pts,colr);
 %plotio_curve(aa,resp(idxx),pts,colr);
end