% Function to find the entropy of the output code space give the input
% code space measure

function [out_ent,in_ent,uni_GC_p2,output_tdex]=trans_ent(in_measure,GC_p1,GC_p2)

in_ent=-sum(in_measure.*log2(in_measure));

% Construct an index of all input codes leading to each output code
% in GC_p2

uni_GC_p2=union(GC_p2,GC_p2);
out_dims=size(uni_GC_p2);
output_tdex=cell(1,out_dims(2));
out_ent=0;

for i=1:out_dims(2)
    
    dex=find(GC_p2==uni_GC_p2(i));
    output_tdex{1,i}=GC_p1(dex);
    out_ent=out_ent-sum(in_measure(dex).*log2(in_measure(dex)));
    
end    

