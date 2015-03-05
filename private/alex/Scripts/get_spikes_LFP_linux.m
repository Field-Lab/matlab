%% get spike data

date='20110908';
cd('/mnt/muench_data/user/alexandra/scripts')
chs=1:60;
chs(15)=[];
for ch=chs
    tic
    ch
    get_spikes_v2(ch,date);
    toc
end


%% get LFP data

date='20121018';
cd('/mnt/muench_data/user/alexandra/scripts')
chs=1:60;
chs(15)=[];
keyWord='HCseq';
for ch=chs
    tic
    ch
    get_LFP_linux(ch,date,keyWord);
    toc
end


%% get LFP data: all, by channel

date='20121207';
cd('/mnt/muench_data/user/alexandra/scripts')
chs=1:60;
for ch=chs
    tic
    ch
    get_LFP_linux_all(ch,date);
    toc
end