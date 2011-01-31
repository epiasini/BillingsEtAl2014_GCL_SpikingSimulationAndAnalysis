% Function to load spike data in a sparse representation

function [spikes]=loadspikes(path,basen,extn,filenos)

reps=size(filenos);
reps=reps(2);
max_spikes=0;

% pre_process: determine number of spikes
% each input file is a list of spike times
% for that cell

for i=1:reps
   
    load([path,basen,num2str(filenos(i)),extn]);
    sp_t=eval([basen,num2str(filenos(i))]);
    
    record_len(i)=max(size(sp_t));
    
    if record_len(i)>max_spikes
        
        max_spikes=record_len(i);
        
    end
    
end

% Use -1 to pad the array: Columns are
% cells and rows are times
spikes=-1*ones(reps,max_spikes);


for i=1:reps
   
    sp_t=eval([basen,num2str(filenos(i))]);
    spikes(i,1:record_len(i))=sp_t;
    
end
 


