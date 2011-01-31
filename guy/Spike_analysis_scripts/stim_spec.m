% Function to find the stimulus specific information
% Guy Billings UCL 2010

function [issi]=stim_spec(clusters,alphabet,pr0,ps,stimuli)

[prgs,prs,pr]=find_dists(pr0,ps,alphabet,stimuli,clusters);
stims=max(stimuli);

% firstly determine the specific information

psgr=(ps'*ones(1,clusters)).*prgs./(pr'*ones(1,stims))';
ker=psgr.*log(psgr);
ker(isnan(ker))=0;
ker=sum(ker);
hsgr=-1*pr.*ker;
hs=-1*ps*log(ps)';
isp=hs-hsgr;

% now the stimulus specific information

for s=1:stims
  issi(s)=prgs(s,:)*isp';
end
    