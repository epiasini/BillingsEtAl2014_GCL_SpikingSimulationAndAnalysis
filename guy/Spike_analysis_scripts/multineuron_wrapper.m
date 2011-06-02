% this is meant to offer the same interface as pdist, so that it could be
% used in the linkage function.

function D = multineuron_wrapper(X)
% X is observations x cells x timepoints
% D is a row vector of length m(m–1)/2, corresponding to pairs of
%  observations in X. The distances are arranged in the order (2,1), (3,1),
%  ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1)).

[obs, cells, timepoints] = size(X);

D = [];

%fprintf(1, '%d/00000', obs*(obs-1)/2);
for h = 1:obs
    for k = h+1:obs
        fprintf(1, '%d/%d\n', size(D, 2), obs*(obs-1)/2);
        D = [D multineuron_distance(squeeze(X(h,:,:)), squeeze(X(k,:,:)))];
    end

%fprintf(1, '\n');
    
end