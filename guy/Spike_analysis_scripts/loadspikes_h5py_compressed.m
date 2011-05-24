function [spikes, stim_number, trial_number] = loadspikes_h5py_compressed(filename, cell_type_index)
    
    hinfo = hdf5info(filename);
    
    stim_number = size(hinfo.GroupHierarchy.Groups, 2);
    trial_number = size(hinfo.GroupHierarchy.Groups(1,1).Groups, 2);
    
    
    
    grc_number = hinfo.GroupHierarchy.Groups(1,1).Groups(1,1).Datasets(1,cell_type_index).Dims(1);
    grc_max_spikes = 0;
    
    for stim=1:stim_number
        for trial=1:trial_number
            stim
            trial
            grc_spike_number = hinfo.GroupHierarchy.Groups(1,stim).Groups(1,trial).Datasets(1,cell_type_index).Dims(2);
           
            if grc_max_spikes < grc_spike_number;       
                grc_max_spikes = grc_spike_number;
            end
        end
    end
    
    spikes = -1 * ones(stim_number*trial_number, grc_number, grc_max_spikes);
    
    for stim=1:stim_number
        for trial=1:trial_number
            grc_spike_number = hinfo.GroupHierarchy.Groups(1,stim).Groups(1,trial).Datasets(1,cell_type_index).Dims(2);
            spikes((stim-1)*trial_number+trial, :, 1:grc_spike_number) = hdf5read(hinfo.GroupHierarchy.Groups(1,stim).Groups(1,trial).Datasets(1,cell_type_index));
        end
    end
    