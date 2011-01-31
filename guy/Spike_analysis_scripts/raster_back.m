% Function to reconstruct a raster plot from the typical code
% associated with a cluster

function [r_data]=raster_back(code,units,s_mode,plot_out)

% unslice the code
dims=size(code);
time_slices=dims(2)/units;
r_data=zeros(units,time_slices);

switch lower(s_mode)
    
    case 'sources'

        for t=1:time_slices
    
            r_data(:,t)=code((t-1)*units+1:t*units);
    
        end    
    
    case't_slice'
        
        for s=1:units
    
            r_data(s,:)=code((s-1)*time_slices+1:s*time_slices);
    
        end  
        
end

if plot_out==1
       figure
       imagesc(r_data');
end 
        
        
