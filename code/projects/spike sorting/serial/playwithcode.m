j=20;
robmean=zeros(1,8);
   if ~isempty(clusters) 
                    like=[];
                    for i=1:size(clusters.means,1)  %the number of clusters
                        if any(clusters.covariances(i,:))   % any is a boolean that sees if there exists a non zero element.  (true if there exists,  falso if doesn't)
                            gg=gmdistribution(clusters.means(i,:),clusters.covariances(i,:));
                            like(:,i)=pdf(gg,prjs');%*clusters.probability(i);  still going electrode by electrode
                        end
                    end

                    if 0 %plot   
                        clf 
                        di=[1 2;1 3;2 3];
                        col=colormap('lines');
                        [junk tlike]=max(like,[],2);
                        for i=1:size(clusters.means,1)
                            for ii=1:3
                                subplot(1,3,ii)
                                hold on
                                tt=find(tlike==i);
                                tt=tt(1:500);
                                plot(prjs(di(ii,1),tt),prjs(di(ii,2),tt),'.','color',col(i,:));
                            end
                        end
                    end
                    
                    [s1 s2]=sort(like,2);  % s2 is the index matrix that keeps track cluster number
                    t_like=s2(:,end);  %identifies the higheset prob cluster for each spike
                    tr_like=s1(:,end)./s1(:,end-1);   % the "gain " in probability  (the length of number of spikes)   should all be greater than 1  bigger the better
                    
                    for i=1:size(clusters.means,1)  %number of clusters
                        t1=find(t_like==i);  %a series of indices

                        cell_id=(j-1)*15+i;  % index cell id
                        
                        index=get_cell_indices(datarun,cell_id);
                        if ~isempty(datarun.spikes{index})
                            datarun=get_contamination(datarun,cell_id);

                            if(datarun.contamination(index)<=thr_contam); % contamination check  brings cell count down... and we get info for each on tcluster
                                robmean(1,i)=robust_mean(log(tr_like(t1)));        %  robust mean of the log of the "gain"
                                 %counting a non-contam cluster
                            end
                        end
                    end
                end