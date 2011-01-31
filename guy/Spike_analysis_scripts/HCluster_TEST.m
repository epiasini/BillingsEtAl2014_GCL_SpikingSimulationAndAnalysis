% Test script to apply hierachical clustering to spike train data
% obtained from multiple trial simulations
% Guy Billings, UCL 2010
%---------------------------------------------------------
% User parameters: (NOTE matlab crashes with only 1 chunk for some reason)

time_win=       250;        % Chunk length (ms)
dt=             2;          % Analysis precision (ms) 
rec_len=        500;       % Observation length (ms)
tau=            20;         % MvR metric tau (ms)
smode=          't_slice';  % slicing mode, time or space
clink=          'ward'; % Clustering linkage method

%---------------------------------------------------------
% Load data:

path=           '/Users/guybillings/project_net_transforms/data/HClust_test/sparse_0.1/'; 
basen_inputs=   'm_spikes_hierarchy_testrep_';
basen_outputs=  'g_spikes_hierarchy_testrep_';
extn=           '_pin0.1_bias-0.0_2010-04-12.14-54-31.dat';
patts=          10;
reps=           10;

observations=patts*reps;

%build file table using a cell array
ti={};
to={};
for p=0:patts-1
    
    for r=0:reps-1
        
        fi={[path,basen_inputs,num2str(r),'_pat',num2str(p),extn]};
        ti=[ti,fi];
        fo={[path,basen_outputs,num2str(r),'_pat',num2str(p),extn]};
        to=[to,fo];

    end
    
end    

spikes_in=loadspikes_compiled(observations,ti);
spikes_ou=loadspikes_compiled(observations,to);

%---------------------------------------------------------

cells_in=size(squeeze(spikes_in(1,:,:)));
cells_in=cells_in(1);
cells_ou=size(squeeze(spikes_ou(1,:,:)));
cells_ou=cells_ou(1);

chunks=floor(rec_len/time_win);
chunked_data_in=zeros(observations,cells_in,chunks,ceil(time_win/dt));
chunked_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

conv_data_in=zeros(observations,cells_in,chunks,ceil(time_win/dt));
conv_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

code_vector_in=zeros(observations,chunks,cells_in*ceil(time_win/dt));
code_vector_ou=zeros(observations,chunks,cells_ou*ceil(time_win/dt));

for o=1:observations
    
    % Chunk the data: Represent spike trains as binary stings
    
    chunked_data_in(o,:,:,:)=chunk(squeeze(spikes_in(o,:,:)),cells_in,time_win,dt,rec_len);
    chunked_data_ou(o,:,:,:)=chunk(squeeze(spikes_ou(o,:,:)),cells_ou,time_win,dt,rec_len);
    
    % Apply MvR metric: Convolve with exponential function
    % NOTE that we could look at clustering of the whole network
    % response and we can look at clustering of individual cell respones
    % (i.e. we dont concatentate the convolved spike trains)
   
    conv_data_in(o,:,:,:)=exp_conv(squeeze(chunked_data_in(o,:,:,:)),dt,tau);
    conv_data_ou(o,:,:,:)=exp_conv(squeeze(chunked_data_ou(o,:,:,:)),dt,tau);
    
    % Slice the data into a single vector for creating of distance matrix
    % matrix having 'observations' rows and 'sources'*'ceil(time_win/dt)'
    % columns
    
    code_vector_in(o,:,:)=slice(squeeze(conv_data_in(o,:,:,:)),smode);
    code_vector_ou(o,:,:)=slice(squeeze(conv_data_ou(o,:,:,:)),smode);
    
end

% Perform hierachical clustering on
% distance matrix with MvR metric

in_dist=zeros(chunks,observations*(observations-1)/2);
ou_dist=zeros(chunks,observations*(observations-1)/2);

in_tree=zeros(chunks,observations-1,3);
ou_tree=zeros(chunks,observations-1,3);

for chunk=1:chunks
    
    % Evaluate distance between code vectors
    
    in_dist(chunk,:)=pdist(squeeze(code_vector_in(:,chunk,:)));
    ou_dist(chunk,:)=pdist(squeeze(code_vector_ou(:,chunk,:)));
    
    % Generate cluster information
    
    in_tree(chunk,:,:)=linkage(squeeze(in_dist(chunk,:)),clink);
    ou_tree(chunk,:,:)=linkage(squeeze(ou_dist(chunk,:)),clink);
    
end    

% Generate array mapping inputs into stimuli classes

stimuli=stim_index(reps,observations);

% Determine information transmitted as a function of the number of clusters
% Compute both the mutual information and the stimulus specific information
% (for a single slice here - could also find average over training slices.
% This quantifies how well the decoder works on the training data

data_tree=in_tree;
chunked_data=chunked_data_in;
data=conv_data_in;
vector=code_vector_in;

train_chunk=1;
ps=ones(1,patts)*1/patts;
[mi,issi]=info_heiracy(observations,ps,stimuli,squeeze(data_tree(train_chunk,:,:)));

% Loop over numbers of clusters and find a decoder in each case using
% train_s slices of data for training.

test_chunk=2;
mi_dec=zeros(1,observations);
issi_dec=zeros(observations,patts);

for opt_clust=2:observations

    o_clustering=cluster(squeeze(data_tree(test_chunk,:,:)),'maxclust',opt_clust);%,'depth',opt_clust);%-(opt_clust-1));
    [conv_book,spike_book]=symbols(opt_clust,o_clustering,...
        squeeze(chunked_data(:,:,test_chunk,:)),squeeze(data(:,:,test_chunk,:)));
    % slice, treating clusters like chunks
    book_vector=slice(conv_book,smode);
    
    % Recalculate mean information retrieval using the unseen data

    [mi_dec(opt_clust),issi_dec(opt_clust,:)]=info_decode(observations,ps,stimuli,book_vector,squeeze(vector(:,test_chunk,:)));
    
end    

 



