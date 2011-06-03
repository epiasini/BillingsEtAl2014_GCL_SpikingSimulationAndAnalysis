function d = multineuron_distance(p, q)

% single observation of a network: cells x (time_win/dt) array

[cells, timepoints] = size(p); 

c = 1;
c_matrix = ones(cells);
c_matrix(~eye(size(c_matrix))) = c;

delta = p - q;

% calculate matrix of inner products between trace differences:
% E_hk = <delta_h|delta_k>

E = delta*delta';

% weight these inner products by taking into account the information on the
% angles between vectors in "cell index space"
weighted_distances = E.*c_matrix;

d = sqrt(sum(sum(weighted_distances)));

end