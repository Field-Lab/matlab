function res = RejectSpaces(sg)
res = '';
for k = 1:length(sg)
    x = sg(k);
   if(isspace(x) == 0)
       res = [res,x];
   end
end
end