% Function to compute the mutual information between inputs and outputs and
% stimulus specific information as a function of the decoder derived 
% from the heirachical clustering.
% Guy Billings, UCL 2010
% output is a 2 column array where the rows are the levels in the heirachy
% column 1 is the number of clusters at that level, column 2 is the mutual
% information. Also given is an array of the stimulus specific information
% for each stimulus.

function mi=info_heiracy(data_dims,ps,stimuli,code_tree)

% Assume each output is unique
%code_list(chunked_data_ou,smode,slicedex);
pr0=1/data_dims;

% Step through the tree and determine psr, and pr at each stage 

mi=zeros(data_dims(1)-1,2);
%issi=zeros(data_dims(1)-1,max(stimuli));

for clusters=2:data_dims(1)
    
    % alphabet is calculated by the clustering algorithm.
    alphabet=cluster(code_tree,'maxclust',clusters);
    mi(clusters-1,1)=clusters;
    mi(clusters-1,2)=mutual_info(clusters,alphabet,pr0,ps,stimuli);
    %issi(clusters-1,:)=stim_spec(clusters,alphabet,pr0,ps,stimuli);

end




 