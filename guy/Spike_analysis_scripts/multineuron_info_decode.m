% Function to compute the mutual information between inputs and outputs and
% stimulus specific information as a function of the decoder derived 
% from the heirachical clustering.
% Guy Billings, UCL 2010
% output is a 2 column array where the rows are the levels in the heirachy
% column 1 is the number of clusters at that level, column 2 is the mutual
% information. Also given is an array of the stimulus specific information
% for each stimulus.

function [mi,issi]=multineuron_info_decode(observations,ps,stimuli,conv_book,conv_data)

[observations, cells, timepoints] = size(conv_data);

[cells, clusters, timepoints] = size(conv_book);

% Assume each output is unique
% code_list(chunked_data_ou,smode,slicedex);
pr0=1/observations;

% alphabet is now calculated by the decoder, not by the 
% clustering algorithm.
%alphabet=decodeFST(observations,clusters,book_vector,data_vector,bdims(2),0);
alphabet=multineuron_decode(observations,clusters,conv_book,conv_data);
mi=mutual_info(clusters,alphabet,pr0,ps,stimuli);
issi=stim_spec(clusters,alphabet,pr0,ps,stimuli); 






  