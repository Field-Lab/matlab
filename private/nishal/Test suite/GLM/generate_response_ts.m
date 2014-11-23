function response=generate_response_ts(mov_params,cell_params)

nTrials=mov_params.nTrials;
s_use=cell_params{1}.stas(:,:,:,end:-1:1);
%s_use(:,:,:,15:end)=0; % DOUBT .. TODO ? Ask EJ
mov=mov_params.mov;
movie_time=size(mov,4);

gen_signals=zeros(movie_time,3);
for col=1:3;
    col
    st_temp=zeros(size(s_use,2),size(s_use,1),1,size(s_use,4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=s_use(:,:,col,itime)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    s_use_new=st_temp;

Filtlen = size(s_use_new,4);
Filtdim1=size(s_use_new,1);
Filtdim2=size(s_use_new,2);

movie_new_len=movie_time;
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=squeeze(mov(:,:,col,:)); % Append zeros before the movie
sz=max(size(mov2,3)-size(s_use_new,4) + 1, 0);

gen_signals(:,col) = reshape(convn(mov2,squeeze(s_use_new(end:-1:1,end:-1:1,1,:)),'valid'),[sz,1]);
end
    %
gen=sum(gen_signals,2);
    
% gen is linear output
postSpikeFilter=cell_params{1}.postSpikeFilter;
postSpikeFilterLen=length(postSpikeFilter);
binsPerFrame=cell_params{1}.binsPerFrame;
binLen=mov_params.refresh/binsPerFrame;
binsTotal=movie_time*binsPerFrame;

totalGen=zeros(nTrials,length(binsTotal),1);
spksGen=zeros(nTrials,length(binsTotal),1);
for iTrial=1:nTrials
PostSpikeGen = zeros(binsTotal,1)+cell_params{1}.tonicDrive;

mov_frame_number=zeros(length(binsTotal),1);

ibincnt=0;
for iframe=1:movie_time
    for ibin=1:binsPerFrame
    ibincnt=ibincnt+1;
    totalGen(iTrial,ibincnt)=exp(gen(iframe) + PostSpikeGen(ibincnt));
    spksGen(iTrial,ibincnt) = double(rand(1)<totalGen(iTrial,ibincnt)*binLen/1000);
    mov_frame_number(ibincnt)=iframe;
    
    futureLen=min(postSpikeFilterLen,binsTotal-ibincnt+1);
    if(spksGen(iTrial,ibincnt)>0)
        PostSpikeGen(ibincnt:ibincnt+futureLen-1) = PostSpikeGen(ibincnt:ibincnt+futureLen-1)+postSpikeFilter(1:futureLen);
    end
    
    end
end
end
response.spksGen=spksGen;
response.mov_frame_number=mov_frame_number;
response.totalGen=totalGen;
end