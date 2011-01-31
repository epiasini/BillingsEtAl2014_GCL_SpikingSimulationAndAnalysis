% Script to produce prototype code plot with k means clustering
% Guy Billings UCL 2009
%--------------------------------------%

%path=           '../data/kmeans_testing/data/Gran_net/PARIS_data/';
%basen_p1=       'MF_3D_Small_';
%basen_p2=       'GrC_3D_Small_';
path=           '../data/kmeans_testing/IF_network/test_data/';
basen_p1=       'CellsA_';
basen_p2=       'CellsB_';
extn=           '.SPIKE_-50.spike';
filenos_p1=     [0:9];
filenos_p2=     [0:9];
s_mode=         't_slice';

time_win=       100; %(ms)
tau=                  5;   %MvR metric (ms)
nhoods=          3;
alpha=             1.1;
cluster_reps=  7;
plot_principle=nhoods;

%--------------------------------------%

sources_p1=size(filenos_p1);
sources_p1=sources_p1(2);
sources_p2=size(filenos_p2);
sources_p2=sources_p2(2);

% Load Data

[spikes_p1]=loadspikes(path,basen_p1,extn,filenos_p1);
[spikes_p2]=loadspikes(path,basen_p2,extn,filenos_p2);

% Perform clustering of input spikes

[codebook_p1,codeindex_p1,members_p1,code_vector_p1,dt_p1,sumD_p1,D_p1,rdi]=...
    code_inout(time_win,nhoods,cluster_reps,spikes_p1,sources_p1,tau,s_mode);

% Perform clustering of output spikes

[codebook_p2,codeindex_p2,members_p2,code_vector_p2,dt_p2,sumD_p2,D_p2,rdo]=...
    code_inout(time_win,nhoods,cluster_reps,spikes_p2,sources_p2,tau,s_mode);

% Convert spike vectors to code lists

dims_p1=size(code_vector_p1);
gamma(1)=1/dims_p1(2);
dims_p2=size(code_vector_p2);
gamma(2)=1/dims_p2(2);
mgamma=min(gamma);
[GC_p1]=toGC(code_vector_p1,alpha,mgamma);
[GC_p2]=toGC(code_vector_p2,alpha,mgamma);

% Estimate population and temporal sparsity
% (find population sparsity within clustering timewindow)

[s_sparse,t_sparse,n_spi_p1,n_spi_p2]=sparsity_est(spikes_p1,spikes_p2,tau);

% Analyse the entropy of input and output codes

mdims=size(union(GC_p1,GC_p1));
in_measure=ones(1,max(size(GC_p1)))/mdims(2);%Assuming equally likely (non-parametric)
[out_ent,in_ent,uni_GC_p2,output_tdex]=trans_ent(in_measure,GC_p1,GC_p2);

% Analyse code separation

[cd,c,d]=clust_div(members_p1,members_p2,codeindex_p1,codeindex_p2,uni_GC_p2,GC_p1,GC_p2);
[rel_total_size,rel_avg_sep,Dmat_p1,Dmat_p2]=clust_stats(sumD_p1,sumD_p2,codebook_p1,codebook_p2);

% plot composition of principle clusters

if plot_principle>0
    
    mags=sort(members_p2,'descend');
    rank(1)=min(find(members_p2==mags(1)));
    rank(2)=min(find(members_p2==mags(2)));
    
    for cluster_p=1:plot_principle
    
        figure
        dc(cluster_p)=cluster_plot(min(find(members_p2==mags(cluster_p))),GC_p1,GC_p2,codeindex_p1,codeindex_p2);
        
    end 
    
end    


    
