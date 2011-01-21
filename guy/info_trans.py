# Script to attempt transmission of activation patterns through the network in the relevant
# threshold regimes. 
# The patterns can be saved.
# Guy Billings, UCL 2010.
##########################
# Load libraries

from java.io import *
from java.lang import Runtime
from math import *
from java.util import Random

import time
import os
import run_support # Companion library to this script

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
simulatorSeed = 4321

sim_path="/home/ucgbgbi/data/optimal_d/GranCellLayer/"
nCprojectfile="GranCellLayer.ncx"
inmod='Input_2' # Name of input to be modified (must be pre-existing in project)
networkMLfilebase = "savedNetworks/NML_"

mbasen='MF_3D_Small_'
msaven='m_spikes.txt'
gbasen='GrC_3D_Small_'
gsaven='g_spikes.txt'
exten='.SPIKE_ -20.spike'
psaven='input_patterns_'
spikes_out=sim_path+"simulations/" # Output path for spike files
sref='test' # Reference string for the simulations

pattern_duration=1000
toplim=10 # Number of patterns
# Bias current for first four threshold points (separately empirically determined)
bias_points=[0,-0.007,-0.012,-0.015]
startpdex=1
endpdex=2
pinc=0.1

# Set node number
node_num=0 # Number from 0

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
p = MpiConfiguration("7core");
p.getHostList().add(MpiHost("localhost",7, 1));
simConfig.setMpiConf(p)

#pm.doLoadNetworkMLAndGenerate(nmlFile,0)
pm.doGenerate("Test_energy_eff", neuroConstructSeed)
simConfig = myProject.simConfigInfo.getSimConfig('Test_energy_eff')  

rng=Random()

##########################
# Generate code  
  
while pm.isGenerating():
  print "Waiting for the project to be generated..."
  time.sleep(2)
  
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
    patttstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)
    pattern_file=open(sim_path+psaven+"_pin"+str(pin)+"_bias"+str(bias_points[node_num])+"_"+patttstr+'.txt',"a+")
    
    for i in range(0, toplim):
  
      now = time.localtime(time.time())
      tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)
      simRef = sref+"_pat"+str(i)+"_pin"+str(pin)+"_bias"+str(bias_points[node_num])+"_"+tstr
  
      id=str(i)
      cell_str=[]
      # Need to change this that patterns where 0 inputs are active are rejected
      for j in range(1,n_mf):
        if rng.nextFloat()<=pin:
            cell_str.append(j)
          
      pattern_file.write(str(cell_str)+'\n')
      simprop_path=sim_path+'simulations/'+simRef+"/simulator.props"
      myProject.generatedElecInputs.reset();

      for k in range(0,len(cell_str)-1):
    
        instr=inmod#+str(k)
        print cell_str[k]
        myProject.generatedElecInputs.addSingleInput(instr,'RandomSpikeTrain','MF_3D_Small',cell_str[k],0,0,None)
      
      #Set the thresholding current
      bias_input = myProject.elecInputInfo.getStim("Bias")
      bias_input.setAmp(NumberGenerator(bias_points[node_num]))
      myProject.elecInputInfo.updateStim(bias_input)
      bias_input = myProject.elecInputInfo.getStim("Bias")
      print bias_input
      
      # Save to a NetworkML file
      pm.saveNetworkStructureXML(myProject, File(sim_path+networkMLfilebase+simRef+".nml"), 0, 0, simConfig.getName(), "Physiological Units")
      print "Network structure saved to file: "+sim_path+networkMLfilebase+simRef+".nml"

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

       
##########################          
# compile spike results 
      
        print "Compiling spikes..."
        run_support.compile_spikes(spikes_out+simRef+'/',id,mbasen,exten,msaven,0,n_mf)
        run_support.compile_spikes(spikes_out+simRef+'/',id,gbasen,exten,gsaven,0,n_gr)
    pattern_file.close()

	  

  


                                                 



