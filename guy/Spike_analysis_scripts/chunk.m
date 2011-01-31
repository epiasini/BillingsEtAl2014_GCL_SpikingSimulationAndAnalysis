% Function to convert spike times in to segments of binary
% vector data and to remove uneeded bins from those vectors

function [data]=chunk(spikes,sources,time_win,dt,rec_len)

% Set the bin size to the
% smallest observed ISI across the whole data set

dims=size(spikes);
ISIs=spikes(:,2:dims(2))-spikes(:,1:dims(2)-1);

% remove negative ISIs (-1 is used as padding in the
% spike arrays) 

ISIs=ISIs(find(ISIs>0));

% chunks are defined by the time window duration
% and the rebinning. 
% Each chunk provides a single code vector for 
% analysis.

chunks=floor(rec_len/time_win);
data=zeros(sources,chunks,ceil(time_win/dt));
 
for chunk=1:chunks
    
    interval_start=(chunk-1)*time_win;
    interval_end=chunk*time_win;
     
    for source=1:sources
        
        valid_times=spikes(source,find(spikes(source,:)<=interval_end));
        valid_times=valid_times(find(valid_times>interval_start));
        valid_times=valid_times-interval_start;
        
        data(source,chunk,ceil(valid_times/dt))=1;
        
    end
    
end
       