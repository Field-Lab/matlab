function res=ei_get_axon(datarun,cell_id,varargin)
% template matching log not normalized seed, normalized derivative(seed) 
%numberstd=8; include=[]; p=2;
%greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('include', []);% 
    p.addParamValue('plot', 0);% 
    p.addParamValue('threshold', 2.5);% 
    p.addParamValue('template', []);%
    p.addParamValue('seed_electrode', []);%
    p.addParamValue('hull', 0);%
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

index=get_cell_indices(datarun,cell_id);
ei=datarun.ei.eis{index};

        %take negativ derivative of max electrode
        if ~isempty(params.template)
            template=params.template(1:end-1);
        else
            if ~isempty(params.seed_electrode)           
                template=-diff(ei(params.seed_electrode,:));
            else   
                [m mm]=max(max(abs(ei')));
                template=-diff(ei(mm,:));
            end
        end


        mcorr=zeros(512,1);
        for i=1:length(datarun.ei.position)
            %[junk,mcorr(i)]=shift_corr(template(2:80),ei(i,2:80));
            [junk,mcorr(i)]=shift_corr(template(1:end),ei(i,2:end));
        end



        electrode_list=find(mcorr<params.threshold);

        electrode_list=union(electrode_list, params.include);


        if params.plot
            plot_ei(datarun,cell_id,'neg_color',[0 0 0],'pos_color',[0 0 0],'fliplr',1,'flipud',1);
            hold on;
            tem=datarun.ei.position(electrode_list,:);
            plot(tem(:,1), tem(:,2),'r*');
        end




        %test on conection
            test_list=[];
            electrode_list_co=[];
            i=1;
            while length(test_list)~=length(electrode_list)
                t1=electrode_list(find(~ismember(electrode_list,test_list)));

                t2=electrode_list(find(is_ei_neighbor(t1(1),electrode_list,datarun.piece.array_id)));
                test_list=union(test_list, t2);

                
                if length(t2)>length(electrode_list_co)
                    electrode_list_co=t2;
                end
                i=i+1;
            end

            electrode_list_co=union(electrode_list_co, params.include);

         %t1=position(electrode_list_co,1);
         %t2=position(electrode_list_co,2);
         %if isequal(t1,ones(size(t1))*t1(1)) & isequal(t2,ones(size(t2))*t2(1))
         if params.hull
             try
                %convex polygon
                po=datarun.position(electrode_list_co,:);
                hull=po(convhull(po(:,1), po(:,2)),:);

                res=find(inpolygon(datarun.position(:,1),datarun.position(:,2),hull(:,1),hull(:,2)));

                if params.plot
                    plot(hull(:,1), hull(:,2));
                end
             catch
                %res=electrode_list_co;
                res=[];
             end
         else
             res=electrode_list_co;
         end

        
        if params.plot
            tem=datarun.ei.position(res,:);
            plot(tem(:,1), tem(:,2),'go');

            %set(gca,'YLim',[-500 500],'XLim',[-995 995],'YTick',[],'XTick',[],'PlotBoxAspectRatio',[2 1 1]);
            set(gca,'YLim',datarun.ei.array_bounds_y,'XLim',datarun.ei.array_bounds_x,'YTick',[],'XTick',[],'PlotBoxAspectRatio',[2 1 1]);
            
            box on;
        end


        
       
        
function [shift,mmin,tem]=shift_corr(data1,data2,fill)
%data1=temp;
%data2=ei(data.channels(cell_list(1))+210,1:80);

%clf;plot(ei(327,10:end-10)');hold on;plot(ei(301,10-shift:end-10-shift)','r');

if nargin==2
    fill=0;
end

t=max(abs(data1));
if t, data1=data1/t; end

t=max(abs(data2));
if t, data2=data2/t; end

tdata1=ones(1,length(data1)*3)*fill;
tdata1(length(data1)+1:length(data1)*2)=data1;

%tdata1=log(tdata1+500)-log(500);
%data2=log(data2+500)-log(500);

for i=1:length(data1)*2
   
    
    tem(i)=sum((tdata1(i:i+length(data1)-1)-data2).^2);%old
    
    %chi2 chi2pdf(mmin,length(data)-1(fit amplitude)-1(fit time))
    %tem(i)=sum( (tdata1(i:i+length(data1)-1)-data2).^2./data2 );
    
    
    %vision
        %tem1=sum(abs(data2));
        %if ~tem1,tem1=1e-16;end
        %tem(i)=sum( (tdata1(i:i+length(data1)-1)/sum(abs(tdata1(i:i+length(data1)-1))) - data2/tem1).^2 );
end 

%[mmin,shift]=min(tem(2:end-2));
%shift=shift-length(data1)-1;

[mmin,shift]=min(tem(35:95));
shift=shift-5;