% Function to load compiled spikes as generated
% by the compilespikes function of the python
% simulation support scripts
% Takes a cell array file table and loads contents of all
% files in table into a single table of spike data
% Guy Billings UCL, 2010

function [spikes]=loadspikes_compiled(observations,filetable)

% Load files into a tensor
maxi=observations;
% Determine maximum spike file size
msize=0;
isize=zeros(maxi,2);
for i=1:maxi
    try
      isize(i,:)=size(load(filetable{i}));
    catch  
       warning(['File ',filetable{i},' failed to load'])
    end    
    if isize(i,1)>msize(1) 
        msize=isize(i,:);
    end
end
% Replace padding 0's with -1's
spikes=-1*ones(maxi,msize(2),msize(1));

for i=1:maxi
   if isize(i,1)~=0      
      spikes(i,:,1:isize(i,1))=load(filetable{i})';   
   end 
end
    
    
    
     
    

 