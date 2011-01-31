% Function to determine the relative change in to spacing and 
% size of clusters.

function [rel_total_size,rel_avg_sep,Dmat_p1,Dmat_p2]=clust_stats(sumD_p1,sumD_p2,codebook_p1,codebook_p2)

dims_p1=size(sumD_p1);
dims_p2=size(sumD_p2);

% overall cluster growth:
rel_total_size=sum(sumD_p2)./sum(sumD_p1);

% determine change to cluster separation:
% find distances between all of the clusters
% and the average separation.
% Generate distance matrix

Dmat_p1=zeros(dims_p1(1));
Dmat_p2=zeros(dims_p2(1));

        
for i=1:dims_p1(1)
            
    for j=1:dims_p1(1)
                           
        Dmat_p1(i,j)=sum((codebook_p1(i,:)-codebook_p1(j,:)).^2);
                
    end
    
end    

for i=1:dims_p2(1)
            
    for j=1:dims_p2(1)
                           
        Dmat_p2(i,j)=sum((codebook_p2(i,:)-codebook_p2(j,:)).^2);
                
    end
    
end  

rel_avg_sep=sum(sum(Dmat_p2))/sum(sum(Dmat_p1));
