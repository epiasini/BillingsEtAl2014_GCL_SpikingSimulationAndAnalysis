% Concatenate spikes 
% to provide the binary vectors for clustering

function [code_vector]=slice(data,mode)

dims=size(data);
chunks=dims(2);
sources=dims(1);
t_slices=dims(3);

% mode controls the way that the code vectors
% are constructed:
% sources: Concatenates states of units end-to-end
% t_slice: Concatenates spike trains in turn

switch lower(mode)
    
    case 't_slice'  

       for c=1:chunks
    
          slices=[];
    
           for source=1:sources
   
             slice=squeeze(data(source,c,:))';
             slices=[slices slice];
        
           end
    
          code_vector(c,:)=slices;
    
       end
        
    case 'sources'
        
        for c=1:chunks
    
          slices=[];
    
           for t_slice=1:t_slices
   
             slice=squeeze(data(:,c,t_slice))';
             slices=[slices slice];
        
           end
    
          code_vector(c,:)=slices;
    
        end
       
end        
        
        
       
       