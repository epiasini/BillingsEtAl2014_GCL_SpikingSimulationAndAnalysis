% Function to compute the mutual information between inputs and outputs and
% stimulus specific information as a function of the decoder derived 
% from the heirachical clustering.
% Guy Billings, UCL 2010
% output is a 2 column array where the rows are the levels in the heirachy
% column 1 is the number of clusters at that level, column 2 is the mutual
% information. Also given is an array of the stimulus specific information
% for each stimulus.

function [mi,issi]=info_decode(data_dims,ps,stimuli,book_vector,data_vector)

observations=size(data_vector);
observations=observations(1);

bdims=size(book_vector);
clusters=bdims(1);

% Assume each output is unique
% code_list(chunked_data_ou,smode,slicedex);
pr0=1/data_dims;

% alphabet is now calculated by the decoder, not by the 
% clustering algorithm.
%alphabet=decodeFST(observations,clusters,book_vector,data_vector,bdims(2),0);
alphabet=decode(observations,clusters,book_vector,data_vector);
mi=mutual_info(clusters,alphabet,pr0,ps,stimuli);
issi=stim_spec(clusters,alphabet,pr0,ps,stimuli); 






  