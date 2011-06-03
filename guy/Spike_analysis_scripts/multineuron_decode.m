% function to perform decoding of inputs using the
% designated set of clusters.

function [alphabet]=multineuron_decode(observations,clusters,conv_book,conv_data)

dist=zeros(1,clusters);
alphabet=zeros(1,observations);
for observation=1:observations
    for cluster=1:clusters
        dist(cluster) = multineuron_distance(squeeze(conv_data(observation,:,:)), squeeze(conv_book(:,cluster,:)));
    end
    md=find(dist==min(dist));
    num_mins=max(size(md));
    if num_mins>1
        alphabet(observation)=round(rand*(num_mins-1))+1;
    else
        alphabet(observation)=md;
    end        
end

     