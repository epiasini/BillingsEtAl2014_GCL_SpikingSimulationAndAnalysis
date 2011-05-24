% Apply hierachical clustering to spike train data
% obtained from multiple trial simulations and evaluate 
% informativeness of coding based on these clusters.
% Guy Billings, UCL 2010
%---------------------------------------------------------
% User parameters: (NOTE matlab crashes with only 1 chunk for some reason)

time_win=       50;        % Chunk length (ms)
dt=             2;          % Analysis precision (ms) 
rec_len=        300;       % Observation length (ms)
tau=            5;         % MvR metric tau (ms)
smode=          't_slice';  % slicing mode, time or space
clink=          'ward'; % Clustering linkage method

%---------------------------------------------------------
% Load data:

path=           '/Users/guybillings/project_net_transforms/data/HClust_test/sparse_bias/output/'; 
basen_inputs=   'm_spikes_hierarchy_test_rep_';
basen_outputs=  'g_spikes_hierarchy_test_rep_';
extn=           '_pin0.1_bias-0.005_2010-04-12.16-49-55.dat';
%patts=          50;
%reps=           10;

hdf5_filename = sprintf('/home/ucbtepi/code/network/data/10_f.6_b%02d.hdf5',bias)

%observations=patts*reps;

% %build file table using a cell array
% ti={};
% to={};
% for p=0:patts-1
    
%     for r=0:reps-1
        
%         fi={[path,basen_inputs,num2str(r),'_pat',num2str(p),extn]};
%         ti=[ti,fi];
%         fo={[path,basen_outputs,num2str(r),'_pat',num2str(p),extn]};
%         to=[to,fo];

%     end
    
% end    

fprintf('loading spiketimes\n')
[spikes_in, stim_number, trial_number]=loadspikes_h5py_compressed(hdf5_filename, 2);
[spikes_ou, stim_number, trial_number]=loadspikes_h5py_compressed(hdf5_filename, 1);

patts = stim_number;
reps = trial_number;

observations = patts*reps;

%---------------------------------------------------------

cells_in=size(squeeze(spikes_in(1,:,:)), 1);
cells_ou=size(squeeze(spikes_ou(1,:,:)), 1);

fprintf('allocating memory\n')
chunks=floor(rec_len/time_win);
chunked_data_in=zeros(observations,cells_in,chunks,ceil(time_win/dt));
chunked_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

conv_data_in=zeros(observations,cells_in,chunks,ceil(time_win/dt));
conv_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

code_vector_in=zeros(observations,chunks,cells_in*ceil(time_win/dt));
code_vector_ou=zeros(observations,chunks,cells_ou*ceil(time_win/dt));

fprintf('chunking data\n')
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
    fprintf('calculating distances: ')
    in_dist(chunk,:)=pdist(squeeze(code_vector_in(:,chunk,:)));
    ou_dist(chunk,:)=pdist(squeeze(code_vector_ou(:,chunk,:)));
    
    % Generate cluster information
    fprintf('generating clustering information: ')
    in_tree(chunk,:,:)=linkage(squeeze(in_dist(chunk,:)),clink);
    ou_tree(chunk,:,:)=linkage(squeeze(ou_dist(chunk,:)),clink);
    
end    

% Generate array mapping inputs into stimuli classes
fprintf('generating array mapping inputs into stimuli classes\n')
stimuli=stim_index(reps,observations);

% Determine information transmitted as a function of the number of clusters
% Compute both the mutual information and the stimulus specific information
% (for a single slice here - could also find average over training slices.
% This quantifies how well the decoder works on the training data
fprintf('mutual information with clustering\n')
data_tree=in_tree;
chunked_data=chunked_data_in;
data=conv_data_in;
vector=code_vector_in;

train_chunk=3;
ps=ones(1,patts)*1/patts;
[mi,issi]=info_heiracy(observations,ps,stimuli,squeeze(data_tree(train_chunk,:,:)));

% Loop over numbers of clusters and find a decoder in each case using
% train_s slices of data for training.

test_chunk=2;
mi_dec=zeros(1,observations);
issi_dec=zeros(observations,patts);

fprintf('MI with the decoder')
for opt_clust=2:observations

    o_clustering=cluster(squeeze(data_tree(train_chunk,:,:)),'maxclust',opt_clust);
    [conv_book,spike_book]=symbols(opt_clust,o_clustering,...
        squeeze(chunked_data(:,:,test_chunk,:)),squeeze(data(:,:,train_chunk,:)));
    % slice, treating clusters like chunks
    book_vector=slice(conv_book,smode);
    
    % Recalculate mean information retrieval using the unseen data

    [mi_dec(opt_clust),issi_dec(opt_clust,:)]=info_decode(observations,ps,stimuli,book_vector,squeeze(vector(:,test_chunk,:)));
    
end    

save(sprintf('result_b%02d',bias), 'mi', 'mi_dec') 



