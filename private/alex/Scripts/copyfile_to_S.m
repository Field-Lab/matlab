flag=1;
last_copied_mcd=1;
last_copied_heka=1;

mcd_orig_path='\\Cin-tm-1ephys\Alexandra\20130302_2\MEA_mcd\';
heka_orig_path='\\Cin-tm-1ephys\Alexandra\20130302_2\HEKA\';

mcd_dest_path='S:\data\alexandra\MEA_data\20130302_2\MEA_mcd\';
heka_dest_path='S:\data\alexandra\MEA_data\20130302_2\HEKA\';
if ~exist(mcd_dest_path,'file')
    mkdir(mcd_dest_path);
end
if ~exist(heka_dest_path,'file')
    mkdir(heka_dest_path);
end


while flag==1
    mcd=dir([mcd_orig_path,'*.mcd']);
    heka=dir([heka_orig_path,'*.phys']);
    a=length(mcd);
    b=length(heka);
    while last_copied_heka<b
        copyfile([heka_orig_path,heka(last_copied_heka).name],[heka_dest_path,heka(last_copied_heka).name])
        last_copied_heka=last_copied_heka+1;
    end
    while last_copied_mcd<a
        copyfile([mcd_orig_path,mcd(last_copied_mcd).name],[mcd_dest_path,mcd(last_copied_mcd).name])
        last_copied_mcd=last_copied_mcd+1;
    end
    k=clock;
    fprintf('Next Step\tHEKA=%d\tMCD=%d\t\t%d:%d\n',b,a,k(4),k(5))
    pause(60)
end