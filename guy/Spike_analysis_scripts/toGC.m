% Function to turn a binary vectors into generalised codes

function [GC]=toGC(data,alpha,gamma)

% data should be in the format of columns for variables and
% rows for observations

dims=size(data);
for i=1:dims(1)
    
    GC(i)=sum(alpha.^(find(data(i,:)==1)*gamma));
    
end