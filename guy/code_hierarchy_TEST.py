# Script to asses how well random patterns spatial patterns are transmitted by the
# network (data is analysed in matlab using hierarchical clustering)
# Guy Billings, UCL 2010.
# -----
#
# Generate a network
# Loop over sparsity
#   Choose random binary dictionary
#   Loop over threshold
#     Loop over pattern repetitions
#       Compile spikes
#       Save spike files, patterns and network in separate dir
#       Clean up
#  
##########################
# Load libraries

from java.io import *
from java.lang import Runtime
from math import *
from java.util import Random

import time
import os
import run_support # Companion library to this script
import random
import subprocess

from ucl.physiol.neuroconstruct.cell.examples import *
from ucl.physiol.neuroconstruct.cell.utils import *
from ucl.physiol.neuroconstruct.project import *
from ucl.physiol.neuroconstruct.cell.compartmentalisation import *
from ucl.physiol.neuroconstruct.project.cellchoice import *
from ucl.physiol.neuroconstruct.neuron import *
from ucl.physiol.neuroconstruct.nmodleditor.processes import *
from ucl.physiol.neuroconstruct.neuroml import *
from ucl.physiol.neuroconstruct.neuroml.hdf5 import *
from ucl.physiol.neuroconstruct.hpc.mpi import *
from ucl.physiol.neuroconstruct.utils import NumberGenerator

##########################
# User defined data

neuroConstructSeed = 1234
simulatorSeed = random.randint(1, 9999)
cores=7

sim_path="/home/ucgbgbi/data/optimal_d/GranCellLayer/"
nCprojectfile="GranCellLayer.ncx"
inmod='Input_2' # Name of input to be modified (must be pre-existing in project)
networkMLfilebase = "NML_"

mbasen='MF_3D_Small_'
msaven='m_spikes_'
gbasen='GrC_3D_Small_'
gsaven='g_spikes_'
exten='.SPIKE_ -20.spike'
psaven='input_patterns_'
spikes_out="/home/ucgbgbi/data/optimal_d/output/" # Output path for spike files and networkML file
sref='hierarchy_test_' # Reference string for the simulations

pattern_duration=400.0
toplim=50 # Number of patterns
repeats=15  # Number of repeat repetitions
# Bias current for first four threshold points (separately empirically determined)
bias_points=[0,5]#-0.02] # Start and finish bias current in units of the increment
binc=-0.005 # Bias increment
startpdex=1
endpdex=9
pinc=0.1

# Set node number
#node_num=0 # Number from 0

##########################
# Initialise from files

projFile = File(sim_path+nCprojectfile)

print "Loading project from file: " + projFile.getAbsolutePath()+", exists: "+ str(projFile.exists())

pm = ProjectManager(None, None)
myProject = pm.loadProject(projFile)

print "Loaded project: " + myProject.getProjectName() 

simConfig = myProject.simConfigInfo.getSimConfig('Test_energy_eff')  
simConfig.setSimDuration(pattern_duration)

# Set up the simulation configuration for parallel exec: ToDo: Write companion function to do
# this and put in run support
p = MpiConfiguration("multicore");
p.getHostList().add(MpiHost("localhost",cores, 1));
simConfig.setMpiConf(p)

#pm.doLoadNetworkMLAndGenerate(nmlFile,0)
pm.doGenerate("Test_energy_eff", neuroConstructSeed)
simConfig = myProject.simConfigInfo.getSimConfig('Test_energy_eff')  

rng=Random()

##########################
# Generate network
  
while pm.isGenerating():
  print "Waiting for the project to be generated..."
  time.sleep(2)

now = time.localtime(time.time())
tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)  
  
# Save to a NetworkML file
nmlFile=File(spikes_out+networkMLfilebase+tstr+".nml")
pm.saveNetworkStructureXML(myProject, File(spikes_out+networkMLfilebase+tstr+".nml"), 0, 0, simConfig.getName(), "Physiological Units")
print "Network structure saved to file: "+spikes_out+networkMLfilebase+tstr+".nml"

n_mf=myProject.generatedCellPositions.getNumberInCellGroup('MF_3D_Small')
n_gr=myProject.generatedCellPositions.getNumberInCellGroup('GrC_3D_Small')
numGenerated = myProject.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(numGenerated)
myProject.neuronSettings.setNoConsole()

if numGenerated > 0:

  print "Generating NEURON scripts..."
  myProject.neuronSettings.setCopySimFiles(1)

##########################  
# Set up and run individual simulations with input generated patterns
# [NOTE: 'run_support.py' should contain a function that patterns the input activation rather than
#  it being done here]
  
for pdex in range(startpdex,endpdex+1):  
    pin=pinc*pdex
    
    now = time.localtime(time.time())
    tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)
    pattern_file=open(spikes_out+psaven+"_pin"+str(pin)+"_"+tstr+'.txt',"a+")
    
    for i in range(0, toplim):

      cell_str=[]
      
      while len(cell_str)==0:
        for j in range(1,n_mf):
          if rng.nextFloat()<=pin:
              cell_str.append(j)
          
      pattern_file.write(str(cell_str)+'\n')
      myProject.generatedElecInputs.reset();

      for k in range(0,len(cell_str)-1):
    
        instr=inmod#+str(k)
        print cell_str[k]
        myProject.generatedElecInputs.addSingleInput(instr,'RandomSpikeTrain','MF_3D_Small',cell_str[k],0,0,None)
      
      for bias in range(bias_points[0],bias_points[1]):
        
        #Set the thresholding current
        bias_input = myProject.elecInputInfo.getStim("Bias")
        bias_input.setAmp(NumberGenerator(bias*binc)) 
        # ToDo: should also set clamp duration to be equal to sim duration
        myProject.elecInputInfo.updateStim(bias_input)
        bias_input = myProject.elecInputInfo.getStim("Bias")
        print bias_input
        
        # Need to regenerate from saved network structure for bias change to take effect
        pm.doLoadNetworkMLAndGenerate(nmlFile,0)
        
        for rep in range(repeats):
        
          simRef = sref+"rep_"+str(rep)+"_pat"+str(i)+"_pin"+str(pin)+"_bias"+str(bias*binc)+"_"+tstr
          simprop_path=sim_path+'simulations/'+simRef+"/simulator.props"
          simulatorSeed = random.randint(1, 999)
          
          myProject.simulationParameters.setReference(simRef)
          myProject.neuronFileManager.generateTheNeuronFiles(simConfig, None, NeuronFileManager.RUN_HOC,simulatorSeed)
          compileProcess = ProcessManager(myProject.neuronFileManager.getMainHocFile())
          compileSuccess = compileProcess.compileFileWithNeuron(0,0)
        
          if compileSuccess:
                       
            pm.doRunNeuron(simConfig)
            tfile = File(simprop_path)
            cnt=1
            print "Simulation running..."
            while tfile.exists() == 0:
              time.sleep(2)
            
            mf_spikes=0
            gr_spikes=0
            # compile spike results 
            print "Compiling spikes..."
            mf_spikes=run_support.compile_spikes(sim_path+'simulations/'+simRef+'/',spikes_out,simRef,mbasen,exten,msaven,0,n_mf)
            gr_spikes=run_support.compile_spikes(sim_path+'simulations/'+simRef+'/',spikes_out,simRef,gbasen,exten,gsaven,0,n_gr)
            
            print 'Active mossy fibre rate: '+str((mf_spikes/(len(cell_str)*pattern_duration))*1000.0)
            print 'Average granule cell rate: '+str((gr_spikes/(n_gr*pattern_duration))*1000.0)
            
            print 'Removing simulation folder.'
            subprocess.call('rm -rf '+sim_path+'/simulations/'+simRef,shell=True)
            
    pattern_file.close()

	  

  


                                                 



