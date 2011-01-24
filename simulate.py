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
from java.lang import System
#from math import *
from java.util import Random

import time
import os
#import run_support # Companion library to this script
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
cores=2

sim_path="/home/eugenio/phd/nC_projects/if_network/"
nCprojectfile="if_network.ncx"
inmod='MF_stim' # Name of input to be modified (must be pre-existing in project)
networkMLfilebase = "NML_"

mbasen='MFs_'
msaven='m_spikes_'
gbasen='GrCs_'
gsaven='g_spikes_'
psaven='input_patterns_'
spikes_out="/home/eugenio/phd/code/network/data/" # Output path for spike files and networkML file
sref='hierarchy_test_' # Reference string for the simulations

pattern_duration=50.0
toplim=2 # Number of patterns
repeats=1  # Number of repeat repetitions

startpdex=1
endpdex=1
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

simConfig = myProject.simConfigInfo.getSimConfig('Default Simulation Configuration')  
simConfig.setSimDuration(pattern_duration)

#pm.doLoadNetworkMLAndGenerate(nmlFile,0)
pm.doGenerate("Default Simulation Configuration", neuroConstructSeed)
simConfig = myProject.simConfigInfo.getSimConfig('Default Simulation Configuration')

rng=Random()

##########################
# Generate network
  
while pm.isGenerating():
    print "Waiting for the project to be generated..."
    time.sleep(2)

now = time.localtime(time.time())
tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)  
  
# Save to a NetworkML file
#nmlFile=File(spikes_out+networkMLfilebase+tstr+".nml")
#pm.saveNetworkStructureXML(myProject, File(spikes_out+networkMLfilebase+tstr+".nml"), 0, 0, simConfig.getName(), "Physiological Units")
#print "Network structure saved to file: "+spikes_out+networkMLfilebase+tstr+".nml"

n_mf=myProject.generatedCellPositions.getNumberInCellGroup('MFs')
n_gr=myProject.generatedCellPositions.getNumberInCellGroup('GrCs')
numGenerated = myProject.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(numGenerated)
myProject.neuronSettings.setNoConsole()

if numGenerated > 0:
    print "Generating NEURON scripts..."
    myProject.neuronSettings.setCopySimFiles(1)

    ##########################  
    # Set up and run individual simulations with input generated patterns
  
    for pdex in range(startpdex,endpdex+1):  
        pin=pinc*pdex
        
        now = time.localtime(time.time())
        tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)
        pattern_file=open(spikes_out+psaven+"_pin"+str(pin)+"_"+tstr+'.txt',"a+")
        
        for i in range(0, toplim):
            myProject.generatedElecInputs.reset()
            cell_str=[]
            
            while len(cell_str)==0:
                for j in range(1,n_mf):
                    if rng.nextFloat()<=pin:
                        cell_str.append(j)
          
            pattern_file.write(str(cell_str)+'\n')
            print "Total inputs "+str(myProject.generatedElecInputs.getNumberSingleInputs(inmod))
            for k in range(0,len(cell_str)):
                print 'adding MF inputs, MF n.' + str(cell_str[k])
                myProject.generatedElecInputs.addSingleInput(inmod,'RandomSpikeTrain','MFs',cell_str[k],0,0,None)  
            print "Total inputs "+str(myProject.generatedElecInputs.getNumberSingleInputs(inmod))
            simRef = tstr + "_" + str(i)
            myProject.simulationParameters.setReference(simRef)
            myProject.neuronFileManager.generateTheNeuronFiles(simConfig, None, NeuronFileManager.RUN_HOC,simulatorSeed)
            compileProcess = ProcessManager(myProject.neuronFileManager.getMainHocFile())
            compileSuccess = compileProcess.compileFileWithNeuron(0,0)
        
            if compileSuccess:
                pm.doRunNeuron(simConfig)
                propsfile_path=sim_path+"simulations/"+simRef+"/simulation.props"
                tfile = File(propsfile_path)
                cnt=1
                print "simref " + str(simRef)
                print propsfile_path
                print "Simulation running..."
                while tfile.exists() == 0:
                    time.sleep(2)
    pattern_file.close()
System.exit(0)
	  

  


                                                 



