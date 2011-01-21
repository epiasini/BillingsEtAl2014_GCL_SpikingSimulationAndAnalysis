# File containing definitions used in order to support the use of neuroconstruct 
# simualtions for batch mode scientific runs.
# Author: Guy Billings, UCL 2010

# Rapid generation of 2d array
def gen_2darray(m, n, initial=0):

    arr = [None] * m
    for i in range(m):
      arr[i] = [initial] * n
    return arr

# Function to take spike locations (i.e. the index of the neuron spiking
# and the time of spiking from memory and write to disk
def compile_spikes(simpath,outpath,id,basen,exten,saven,start_num,end_num):
    
    max_spikes=0
    tot_spikes=0
    for file_num in range(start_num,end_num):
      file=simpath+basen+str(file_num)+exten
      print file
      fh = open(file,"r")
      spikes=fh.read()
      num_spikes=spikes.count("\n")
      tot_spikes=tot_spikes+num_spikes
      if num_spikes>max_spikes:
        max_spikes=num_spikes  
      fh.close()
    print str(tot_spikes/(end_num-start_num))+' spikes generated on average.'
      
    st=gen_2darray(max_spikes,(end_num-start_num))
    for file_num in range(start_num,end_num): 
      file=simpath+basen+str(file_num)+exten
      fh = open(file,"r")
      line_no=0
      for line in fh:
        if line!='\n':
          st[line_no][file_num]=float(line)
        line_no=line_no+1
      fh.close() 
      
    file=outpath+saven+id+".dat"
    fo=open(file,"w")    
    out_str=''
    for outline in range(1,max_spikes):
      for file in range(0,(end_num-start_num)):
        out_str=out_str+str(st[outline][file])+' '
      out_str=out_str+'\n'
    fo.write(out_str)
    fo.close()  