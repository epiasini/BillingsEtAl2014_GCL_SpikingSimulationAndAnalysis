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

working_dir = sprintf('/home/ucbtepi/code/network/data/f.5_20_-20/s%.2f', scale);
hdf5_filename = sprintf('%s/20_f.5_s%.2f_b%02d.hdf5', working_dir, scale, bias)

fprintf('loading spiketimes\n')
[spikes_in, stim_number, trial_number]=loadspikes_h5py_compressed(hdf5_filename, 2);
[spikes_ou, stim_number, trial_number]=loadspikes_h5py_compressed(hdf5_filename, 1);

patts = stim_number;
reps = trial_number;

observations = patts*reps;

%---------------------------------------------------------

cells_ou=size(squeeze(spikes_ou(1,:,:)), 1);

fprintf('allocating memory\n')
chunks=floor(rec_len/time_win);
chunked_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

conv_data_ou=zeros(observations,cells_ou,chunks,ceil(time_win/dt));

fprintf('chunking data\n')
for o=1:observations
    
    % Chunk the data: Represent spike trains as binary stings
    
    chunked_data_ou(o,:,:,:)=chunk(squeeze(spikes_ou(o,:,:)),cells_ou,time_win,dt,rec_len);
    
    % Apply MvR metric: Convolve with exponential function
    % NOTE that we could look at clustering of the whole network
    % response and we can look at clustering of individual cell respones
    % (i.e. we dont concatentate the convolved spike trains)
   
    conv_data_ou(o,:,:,:)=exp_conv(squeeze(chunked_data_ou(o,:,:,:)),dt,tau);
    
end

% Perform hierachical clustering on
% distance matrix with MvR metric

ou_dist=zeros(chunks,observations*(observations-1)/2);

ou_tree=zeros(chunks,observations-1,3);

for chunk=1:chunks
    
    % Evaluate distance between output signals
    fprintf('calculating distances: ')
    ou_dist(chunk,:) = multineuron_wrapper(squeeze(conv_data_ou(:,:,chunk,:)));

    % Generate cluster information
    fprintf('generating clustering information: ')
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
data_tree=ou_tree;
chunked_data=chunked_data_ou;
data=conv_data_ou;

train_chunk=3;
ps=ones(1,patts)*1/patts;
[mi,issi]=info_heiracy(observations,ps,stimuli,squeeze(data_tree(train_chunk,:,:)));

% Loop over numbers of clusters and find a decoder in each case using
% train_s slices of data for training.

test_chunk=2;
mi_dec=zeros(1,observations);
issi_dec=zeros(observations,patts);

fprintf('MI with the decoder')
% prellocate space for distance matrix and its mask. Initialise mask
distance_matrix = zeros(observations, observations+observations-1);
relevant_columns = zeros(1, observations+observations-1);

for opt_clust=observations:-1:1
    opt_clust
    % create codebook by calculating cluster centers
    o_clustering=cluster(squeeze(data_tree(train_chunk,:,:)),'maxclust',opt_clust);
    [conv_book,spike_book]=symbols(opt_clust,o_clustering,...
        squeeze(chunked_data(:,:,train_chunk,:)),squeeze(data(:,:,train_chunk,:)));
    
    % prepare matrix of distances between observations and codewords
    clust_index = observations+(observations - opt_clust);
    if opt_clust==observations
        for observation=1:observations
            for clust=1:opt_clust
                distance_matrix(observation,clust) =  multineuron_distance(squeeze(data(observation,:,test_chunk,:)), squeeze(conv_book(:,clust,:)));
            end
        end
        relevant_columns(1:observations) = 1;
    else
        for observation=1:observations
            distance_matrix(observation,clust_index) =  multineuron_distance(squeeze(data(observation,:,test_chunk,:)), squeeze(conv_book(:,opt_clust,:)));
        end
        joined = squeeze(data_tree(train_chunk,observations-opt_clust,1:2));
        relevant_columns(joined) = 0; % ignore clusters that don't exist anymore
        relevant_columns(clust_index) = 1; % take into account the new cluster
    end
    
    % Recalculate mean information retrieval using the unseen data
    [mi_dec(opt_clust),issi_dec(opt_clust,:)]=multineuron_info_decode(observations,ps,stimuli,conv_book,squeeze(chunked_data(:,:,test_chunk,:)), distance_matrix(:,logical(relevant_columns)));
    
end

out_filename = sprintf('%s/result_b%02d', working_dir, bias)
save(out_filename, 'mi', 'mi_dec') 



