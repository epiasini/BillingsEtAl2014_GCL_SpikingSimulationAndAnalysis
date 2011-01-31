function stimuli=stim_index(reps,observations)

stimuli=zeros(1,observations);
classes=observations/reps;

for class=1:classes
    stimuli(reps*(class-1)+1:(class)*reps)=class;
end
    




