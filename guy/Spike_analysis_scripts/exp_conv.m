% Function to convolve binary spike train chunks with an exponential
% for application of the van Rossum spike train distance
% Guy Billings 2010, UCL

function rdata2=exp_conv(rdata,dt,tau)

rsize=size(rdata);
sources=rsize(1);
chunks=rsize(2);
bins=rsize(3);
rdata2=zeros(rsize);  

for chunk=1:chunks
        
       for source=1:sources
           
           spikes=find(rdata(source,chunk,:)==1);
       
           if isempty(spikes)==0
               
             for spike=1:max(size(spikes))  
       
               single_spike=zeros(1,1,bins);
               kern=exp(-[0:(bins-spikes(spike))]/(tau/dt));
               single_spike(1,1,spikes(spike):bins)=kern;
           
               rdata2(source,chunk,:)=rdata2(source,chunk,:)+single_spike;
               
             end  
           
           end
           
       end   
       
       % Take account of the across-chunk residual values
       if chunk~=chunks
           residual=squeeze(rdata2(:,chunk,bins))*exp(-[0:bins-1]/(tau/dt));
           rdims=size(residual);
           indat=zeros(rdims(1),1,rdims(2));
           indat(:,1,:)=residual;
           rdata2(:,chunk+1,:)=indat;
       end    
         
            
       
       
end   
    