% Function to determine the temporal sparsity (frequency) and the 
% population sparsity. Temporal sparsity is found in blocks that are 
% the same duration as the data chunks. Population sparsity is found
% within the clustering timewindow.

function [s_sparse,t_sparse,n_spi_p1,n_spi_p2]=sparsity_est(spikes_p1,spikes_p2,tau)

sdim1=size(spikes_p1);
sdim2=size(spikes_p2);

max_time_p1=max(max(spikes_p1));
max_time_p2=max(max(spikes_p2));

t_sparse(1)=max(size(spikes_p1))/(max_time_p1/1000);
t_sparse(2)=max(size(spikes_p2))/(max_time_p2/1000);

windows(1)=ceil(max_time_p1/tau);
windows(2)=ceil(max_time_p2/tau);

n_spi_p1=zeros(1,windows(1));
n_spi_p2=zeros(1,windows(2));
%find spikes in the timewindow and 
%count each neuron that spikes
%Count a single neuron a mamximum of 1 time

for ds=1:2

    for w=1:windows(ds)
        
        time_int_st=(w-1)*tau;
        time_int_en=w*tau;
        
        if ds==1
        
        for unit=1:sdim1(1)
            
            spi=spikes_p1(unit,spikes_p1(unit,:)>=0);
            slist=spi(spi<time_int_en);
            slist=slist(slist>time_int_st);
            if isempty(slist)==0
            
               n_spi_p1(w)=n_spi_p1(w)+1;
               
            end
            
        end
        
        end
        
        if ds==2

        for unit=1:sdim2(1)
            
            spi=spikes_p2(unit,spikes_p2(unit,:)>=0);
            slist=spi(spi<time_int_en);
            slist=slist(slist>time_int_st);
            if isempty(slist)==0
            
               n_spi_p2(w)=n_spi_p2(w)+1;
               
            end
            
        end
        
        end
            
    end
    
end    
   
s_sparse(1)=mean(n_spi_p1);
s_sparse(2)=mean(n_spi_p2);
        
        