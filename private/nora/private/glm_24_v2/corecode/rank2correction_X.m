%%%% Should match Trainpars.lgrad
%%% LOAD UP THE 

clear movieX
[spaceidx1 timeidx1 spaceidx2 timeidx2] = get_sep_filt_idcesAH(0,Basepars);


space_filt1 = p(spaceidx1+1); % nspace x 1
space_filt2 = p(spaceidx2+1); % nspace x 1
% Get the 2 temporal filters for this neuron
temp_filt1 = p(timeidx1+1); % ntime x 1
temp_filt2 =p(timeidx2+1); %



moveoff_front = 101;

%nspace by time in stim frames
smallmovie = Stimpars.movie_ROI(:,1:moveoff_front + 99);
lgrad_part1 = Trainpars.lgrad{1}([spaceidx1 timeidx1 spaceidx2 timeidx2], moveoff_front:moveoff_front + 99  ) ; %% this is kfilter by time in stim frames




X1 = zeros( length(spaceidx1)  + length(timeidx1)  , size(smallmovie,2) );


%{
%%% Convovled with Time for the Space gradient vectors
for i = 1 : size(smallmovie,2)
    timebegin1 = max(1 , i +1  - length(temp_filt1));
    t_filt1 = temp_filt1(1 : (i - timebegin1) + 1 );
    tobeconvolved1 = smallmovie(   :, timebegin1 : i  );
    movieX( (spaceidx1-1) , i ) =    (tobeconvolved1*t_filt1);
    
    timebegin2 = max(1 , i +1  - length(temp_filt2));
    t_filt2 = temp_filt2(1 : (i - timebegin2) + 1 );
    tobeconvolved2 = smallmovie(   :, timebegin2 : i  );
    movieX( (spaceidx2-1) , i ) =    (tobeconvolved2*t_filt2);
end
%}

%%% Convovled with Time for the Space gradient vectors
for i = moveoff_front : moveoff_front + size(lgrad_part1,2)-1
    timebegin = (i-length(temp_filt1) ) + 1;
    tobeconvolved_pre = smallmovie( :, timebegin:(i) );
    tobeconvolved = fliplr(tobeconvolved_pre);
    movieX( spaceidx1 ,i) =    (tobeconvolved*temp_filt1);
    movieX( timeidx1 , i ) =    (tobeconvolved'*space_filt1);
    movieX( spaceidx2 ,i) =    (tobeconvolved*temp_filt2);
    movieX( timeidx2 , i ) =    (tobeconvolved'*space_filt2);
  %  timebegin2 = max(1 , i - length(temp_filt2));
  %  t_filt2 = temp_filt2(1 : (i - timebegin2) );
  %  tobeconvolved2 = smallmovie(   :, timebegin2 : i-1  );
  %  movieX( (spaceidx2 , i ) =    (tobeconvolved2*t_filt2);
end
figure; subplot(2,1,1); plot(Stimpars.dt*movieX(1,moveoff_front:end)); subplot(2,1,2); plot(lgrad_part1(1,:))
figure; subplot(2,1,1); plot(Stimpars.dt*movieX(226,moveoff_front:end)); subplot(2,1,2); plot(lgrad_part1(226,:))
figure; subplot(2,1,1); plot(Stimpars.dt*movieX(260,moveoff_front:end)); subplot(2,1,2); plot(lgrad_part1(260,:))
figure; subplot(2,1,1); plot(Stimpars.dt*movieX(509,moveoff_front:end)); subplot(2,1,2); plot(lgrad_part1(509,:))
max(max(abs( Stimpars.dt * movieX(1:235,moveoff_front:end) -lgrad_part1(1:235,:) ) ))
max(max(abs( Stimpars.dt * movieX(:,moveoff_front:end) -lgrad_part1  )))


%%% Convolved with Space for the Time gradient vector
%{
for i = moveoff_front : moveoff_front + size(lgrad_part1,2)-1
    timebegin = (i-length(temp_filt1) ) + 1;
    tobeconvolved1 = smallmovie( :, timebegin:(i) );
    movieX( timeidx1 , i ) =    (tobeconvolved1'*space_filt1);
    
    timebegin2 = max(1 , i +1  - length(temp_filt2));
    t_filt2 = temp_filt2(1 : (i - timebegin2) + 1 );
    tobeconvolved2 = smallmovie(   :, timebegin2 : i  );
    movieX( (timeidx2(1:size(tobeconvolved2,2)))-1 , i ) =    (tobeconvolved2'*space_filt2);
end
%%%%%%

a=(find(Stimpars.dt *movieX - Trainpars.lgrad{1}(:,1:1000) ));

size(a)
plot(a(1:1000))
%}