% Function to generate clusters in binary code space

function [codebook,codeindex,members,sumD,D,rdata2]=codeQ...
    (input_data,nhoods,reps,method,tau,dt,sources,s_mode)

% use k-means to cluster the data:
% input data must be in format of rows as observations
% and columns as variables

dims=size(input_data);

switch lower(method)
    
    case 'manhat' % standard kmeans with Manhattan distance metric'

    [clusters,codebook,sumD,D]=kmeans...
        (input_data,nhoods,'distance','cityblock','replicates',reps,'emptyaction','drop');
    rdata2=0;
    
    case 'vanros'% kmeans with van Rossum metric
    
    % Filter binary vectors with exponential kernel and cluster with
    % squared Euclidean distance giving van Rossum spike metric.
    
    rsize=[sources,dims(1),dims(2)/sources];
    rdata2=zeros(rsize);  
    
    for chunk=1:dims(1)
        
       % Unslice the input data to apply the kernel
       [rdata]=raster_back(input_data(chunk,:),sources,s_mode,0);
            

       for source=1:sources
           
           spikes=find(rdata(source,:)==1);
       
           if isempty(spikes)==0
               
             for spike=1:max(size(spikes))  
       
               single_spike=zeros(1,1,rsize(3));
               kern=exp(-[0:(rsize(3)-spikes(spike))]/(tau/dt));
               single_spike(1,1,spikes(spike):rsize(3))=kern;
           
               rdata2(source,chunk,:)=rdata2(source,chunk,:)+single_spike;
               
             end  
           
           end
           
       end    
       
    end   
    
    % Re-slice the data 
    input_data2=slice(rdata2,s_mode);
    
    % Perform clustering of the filtered trains
    [clusters,codebook,sumD,D]=kmeans...
        (input_data2,nhoods,'distance','sqEuclidean','replicates',reps,'emptyaction','drop');
    
    D=D/tau;
    sumD=sumD/tau;
    
    %case 'spatiotemp'
    % Use a Parzen window to cluster the data? Choice of the window is
    % non-trivial: Possibly learn it?
    
end    

codeindex=zeros(nhoods,dims(1));

for i=1:nhoods
    contents=find(clusters==i);
    tmp=size(contents);
    members(i)=tmp(1);
    codeindex(i,1:members(i))=contents;
end    



