% Basic function to cluster spike vectors with k-means
% from spike data (array of spike times (see loadspikes.m)

function [codebook,codeindex,members,code_vector,dt,sumD,D,rd]=code_inout...
    (time_win,nhoods,cluster_reps,spikes,sources,tau,s_mode)

% Look at data in chunks of duration of time_win
% Regard each chunk as a 'repeat'

[data,chunks,dt]=chunk(spikes,sources,time_win);

[code_vector]=slice(data,s_mode);

if tau>0
    
    method='vanros';
    
else
    
    method='manhat';
    
end    

[codebook,codeindex,members,sumD,D,rd]=codeQ...
    (code_vector,nhoods,cluster_reps,method,tau,dt,sources,s_mode);



        
        