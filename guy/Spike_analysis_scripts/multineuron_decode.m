% function to perform decoding of inputs using the
% designated set of clusters.

function [alphabet]=multineuron_decode(observations,clusters, distance_matrix)

alphabet=zeros(1,observations);
mindist = min(distance_matrix,[],2);
for observation=1:observations
    distances = distance_matrix(observation,:);
    md=find(distances==mindist(observation));
    num_mins=max(size(md));
    if num_mins>1
        alphabet(observation)=round(rand*(num_mins-1))+1;
    else
        alphabet(observation)=md;
    end        
end

     