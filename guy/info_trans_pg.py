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
##########################
# User defined data

neuroConstructSeed = 1234
simulatorSeed = 4321

nCprojectfile="/home/guy/test_par/GranCellLayer/GranCellLayer.ncx"
nCnetworkfile="/home/guy/test_par/GranCellLayer/savedNetworks/Net_08-Feb-10_15-32-28.nml"
inmod='Input_2' # Name of input to be modified (must be pre-existing in project)

mbasen='MF_3D_Small_'
msaven='m_spikes.txt'
gbasen='GrC_3D_Small_'
gsaven='g_spikes.txt'
exten='.SPIKE_ -20.spike'
spikes_out='/home/guy/test_par/GranCellLayer/simulations/' # Output path for spike files
sref='d_test_' # Reference string for the simulations

toplim=3 # Number of patterns
pin=0.5   # Prob. input is activated in pattern

##########################
# Initialise from files

projFile = File(nCprojectfile)
nmlFile = File(nCnetworkfile)

print "Loading project from file: " + projFile.getAbsolutePath()+", exists: "+ str(projFile.exists())

myProject = Project.loadProject(projFile, Project.getDummyProjectEventListener())
pm = ProjectManager(None, None)
pm.setCurrentProject(myProject)

print "Loaded project: " + myProject.getProjectName() 

simConfig = myProject.simConfigInfo.getSimConfig('Test_energy_eff')  

# Set up the simulation configuration for parallel exec: ToDo: Write companion function to do
# this and put in run support
p = MpiConfiguration("bernal");
p.getHostList().add(MpiHost("localhost",1, 1));
simConfig.setMpiConf(p)

pm.doLoadNetworkMLAndGenerate(nmlFile,0)
n_mf=myProject.generatedCellPositions.getNumberInCellGroup('MF_3D_Small')
n_gr=myProject.generatedCellPositions.getNumberInCellGroup('GrC_3D_Small')

rng=Random()

##########################
# Generate code  
  
while pm.isGenerating():
  print "Waiting for the project to be generated..."
  time.sleep(2)
      
numGenerated = myProject.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(numGenerated)
myProject.neuronSettings.setNoConsole()

quit()
if numGenerated > 0:

  print "Generating NEURON scripts..."
      
  myProject.neuronSettings.setCopySimFiles(1)

##########################  
# Set up and run individual simulations with input generated patterns
# [NOTE: 'run_support.py' should contain a function that patterns the input activation rather than
#  it being done here]

  for i in range(1, toplim):
    id=str(i)
    cell_str=[]

    for j in range(1,n_mf):
      if rng.nextFloat()<=pin:
          cell_str.append(j)
           
    now = time.localtime(time.time())
    tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)
    simRef = sref+str(i)+"_"+tstr
    simprop_path='/home/guy/test_par/GranCellLayer/simulations/'+simRef+"/simulator.props"
    print cell_str
    myProject.generatedElecInputs.reset();
    
    for k in range(0,len(cell_str)-1):
    
      instr=inmod#+str(k)
      print cell_str[k]
      myProject.generatedElecInputs.addSingleInput(instr,'RandomSpikeTrain','MF_3D_Small',cell_str[k],0,0,None)

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
      

	  

  


                                                 



