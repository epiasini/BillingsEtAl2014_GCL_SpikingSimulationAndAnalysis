% script to test decodeFST

load test.mat
debug=0;
clusters=10;
bdims=size(book_vector);
bdims=bdims(2);

alphabetFST=decodeFST(observations,clusters,book_vector,squeeze(vector(:,test_chunk,:)),bdims,debug)