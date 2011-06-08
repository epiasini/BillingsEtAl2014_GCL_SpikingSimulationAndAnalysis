% function to determine the joint and conditional response stimulus
% distributions given a candidate code.

function [prgs,prs,pr]=find_dists(pr0,ps,alphabet,stimuli,clusters)

% Construct the stimulus response joint under the alphabet
stims=max(stimuli);
prgs=zeros(stims,clusters);
prs=zeros(stims,clusters);
for stim=1:stims 
    
    % probability of the reponse given the stimulus

    obs=find(stimuli==stim);
    prgs(stim,:)=hist(alphabet(obs),[1:clusters]);
    prgs(stim,:)=prgs(stim,:)/max(size(obs));
    prs(stim,:)=prgs(stim,:)*ps(stim);
end
pr=zeros(1,clusters);
for cluster=1:clusters
    occurrences = find(alphabet==cluster);
    if ~isempty(occurrences)
        num_r=max(size(occurrences));
        pr(cluster)=num_r*pr0;
    else
        % in this case, size(occurrences) would be [1 0]...
        pr(cluster) = 0;
    end
     
end  