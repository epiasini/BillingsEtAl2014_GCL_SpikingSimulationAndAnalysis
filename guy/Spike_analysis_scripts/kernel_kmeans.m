% Function to impliment kmeans in a more cutomisable format
% such that the clusters are formed in a feature space having a custom
% kernel (see p275, Kernel methods fo pattern recognition, Shaw-Taylot &
% Chrsitiani)
% Clustering using the van Rossum spike train distance metric can thus be
% achieved by means of supplying an exponential kernel.
% Guy Billings UCL 2009

function []=kernel_kmeans()


