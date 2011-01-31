% Function to calculate the mutual information between the stimulus 
% classes and the observed responses given some output alphabet

function mi=mutual_info(clusters,alphabet,pr0,ps,stimuli)

[prgs,prs,pr]=find_dists(pr0,ps,alphabet,stimuli,clusters);
 
dm=ps'*pr; 
mi=prs.*log2(prs./dm);
mi(find(isnan(mi)))=0;
mi=sum(sum(mi)); 