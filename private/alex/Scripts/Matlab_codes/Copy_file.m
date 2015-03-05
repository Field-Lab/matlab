i=1006;
while flag
    tic
%     if check
%         a=dir(mainPath);        
%         display('GOT it, all OK')
%         size(a)
%         heka=dir(hekaPath);        
%         check=0;
%     else   
%         if mod(i,10)==0
%             check=1;
%         end
%         copyfile(strcat(mainPath(1:end-7),'HEKA\',heka(i/2+1).name),strcat(strcat('F:\20110518\HEKA\',heka(i/2+1).name)));        
        copyfile(strcat(mainPath,'\',a(i).name),strcat(strcat('F:\20110518\MEA_mcd\',a(i).name)));        
        i=i+2;
%     end

    toc
end
