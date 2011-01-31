% Function to calculate a measure of the cluster divergence. Gives
% -1 for convergent clusters (i.e. overall, codes are remapped into
% few clusters) or 1 for divergent clusters. The measure is essentially
% the difference of the clustering entropy in the forward direction 
% and the entropy in the backward direction, expressed as a fraction
% of the maximal clustering entropy in each direction.

function [cd,c,d]=clust_div(members_p1,members_p2,codeindex_p1,codeindex_p2,uni_GC_p2,GC_p1,GC_p2)

members_p1=members_p1(find(members_p1~=0));
members_p2=members_p2(find(members_p2~=0));

dims_p1=size(members_p1);
dims_p1=dims_p1(2);
dims_p2=size(members_p2);
dims_p2=dims_p2(2);

% Find the mapping from input codes to codes in the ouput cluster, 
% given the fraction of the inputs for each input cluster that are 
% found in each output cluster.

number_found=zeros(dims_p2,dims_p1);

for oc=1:dims_p2
    
 for code=1:members_p2(oc)
     
     for ic=1:dims_p1 
         
        if isempty(find(GC_p1(codeindex_p1(ic,1:members_p1(ic)))==GC_p1(codeindex_p2(oc,code))))==0
            
            number_found(oc,ic)=number_found(oc,ic)+1;
            
        end
        
     end
     
 end
 
end 

%KL Divergence with the 'ideal coding' 
%[show that uniform dist is optimal for pairwise 
%yes no comparisons?]

zi=members_p1/sum(members_p1);                                       
zo=members_p2/sum(members_p2);  

c=sum(zi.*log(zi/(1/dims_p1)));
d=sum(zo.*log(zo/(1/dims_p2)));
        
cd=d-c;

% Find divergence for each input cluster
%zi=members_p1/sum(members_p1);                                       
%zo=members_p2/sum(members_p2);  
%denomi=log(1/dims_p2);
%denomo=log(1/dims_p1);
%for ic=1:dims_p1
    
%    number_found(find(number_found(:,ic)==0),ic)=sum(members_p1(ic));
%    f=(number_found(:,ic)/members_p1(ic)); 
                                                                  
%    d(ic)=sum(f.*log(f))/denomi;

%    number_found(find(number_found(:,ic)==members_p1(ic)),ic)=0;
    
%end 

%d=sum(zi.*d);

% Find convergence for each output cluster

%for oc=1:dims_p2
    
%    number_found(oc,find(number_found(oc,:)==0))=sum(members_p2(oc));
%    f=(number_found(oc,:)/members_p2(oc))*(members_p2(oc)/sum(members_p2));
    
%    c(oc)=sum(f.*log(f))/denomo;

%    number_found(oc,find(number_found(oc,:)==members_p2(oc)))=0;
%end    

%c=sum(zo.*c);