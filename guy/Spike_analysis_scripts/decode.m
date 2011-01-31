% function to perform decoding of inputs using the
% designated set of clusters.

function [alphabet]=decode(observations,clusters,book_vector,data_vector)

dist=zeros(1,clusters);
alphabet=zeros(1,observations);
for observation=1:observations
    for cluster=1:clusters      
       dist(cluster)=sqrt(sum((book_vector(cluster,:)-data_vector(observation,:)).^2));
    end
    md=find(dist==min(dist));
    num_mins=max(size(md));
    if num_mins>1
        alphabet(observation)=round(rand*(num_mins-1))+1;
    else
        alphabet(observation)=md;
    end        
end

     