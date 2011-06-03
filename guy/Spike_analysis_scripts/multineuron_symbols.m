% function to find the symbols associated with a candidate coding, 
% i.e. the cluster averages for each candidtate cluster
% Find symbols both in feature space and in spike train space.

function [conv_book,spike_book]=multineuron_symbols(clusters,clustering,chunked_data,conv_data)

[observations, cells, timepoints]=size(chunked_data);
conv_book = zeros(clusters, cells, timepoints);
conv_book=zeros(data_dims(2),clusters,data_dims(3));
spike_book=zeros(data_dims(2),clusters,data_dims(3));

for cluster=1:clusters
    
    obs_in_cluster=find(clustering==cluster);
    if max(size(obs_in_cluster))>1
        convadd=reshape(sum(conv_data(obs_in_cluster,:,:)),data_dims(2),1,data_dims(3));
        spikeadd=reshape(sum(chunked_data(obs_in_cluster,:,:)),data_dims(2),1,data_dims(3));
    else 
        convadd=reshape(conv_data(obs_in_cluster,:,:),data_dims(2),1,data_dims(3));
        spikeadd=reshape(chunked_data(obs_in_cluster,:,:),data_dims(2),1,data_dims(3));
    end 
    conv_book(:,cluster,:)=conv_book(:,cluster,:)+convadd;
    conv_book(:,cluster,:)=conv_book(:,cluster,:)/max(size(obs_in_cluster));
    spike_book(:,cluster,:)=spike_book(:,cluster,:)+spikeadd;
    spike_book(:,cluster,:)=spike_book(:,cluster,:)/max(size(obs_in_cluster));
    
end    
 
  
 