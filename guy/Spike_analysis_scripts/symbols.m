% function to find the symbols associated with a candidate coding, 
% i.e. the cluster averages for each candidtate cluster
% Find symbols both in feature space and in spike train space.

function [conv_book,spike_book]=symbols(relevant_codewords,clustering,chunked_data,conv_data)
% TODO: this could be optimized by keeping a record of the cluster centers that have already been calculated at previous steps of the hierarchical procedure. 
[obs cells timepoints] = size(chunked_data); % obs cells timepoints
conv_book=zeros(cells,obs+obs-1,timepoints);
spike_book=zeros(cells,obs+obs-1,timepoints);

for cluster=find(relevant_codewords)
    
    obs_in_cluster=find(clustering==cluster);
    if max(size(obs_in_cluster))>1
        convadd=reshape(sum(conv_data(obs_in_cluster,:,:)),cells,1,timepoints);
        spikeadd=reshape(sum(chunked_data(obs_in_cluster,:,:)),cells,1,timepoints);
    else 
        convadd=reshape(conv_data(obs_in_cluster,:,:),cells,1,timepoints);
        spikeadd=reshape(chunked_data(obs_in_cluster,:,:),cells,1,timepoints);
    end 
    conv_book(:,cluster,:)=conv_book(:,cluster,:)+convadd;
    conv_book(:,cluster,:)=conv_book(:,cluster,:)/max(size(obs_in_cluster));
    spike_book(:,cluster,:)=spike_book(:,cluster,:)+spikeadd;
    spike_book(:,cluster,:)=spike_book(:,cluster,:)/max(size(obs_in_cluster));
    
end    
 
  
 