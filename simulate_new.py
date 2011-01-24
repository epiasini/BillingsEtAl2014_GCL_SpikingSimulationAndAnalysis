# -*- coding: utf-8 -*-
import random
import time
import os.path

from java.lang import System
from java.io import File

from ucl.physiol.neuroconstruct.project import ProjectManager
from ucl.physiol.neuroconstruct.utils import NumberGenerator
from ucl.physiol.neuroconstruct.nmodleditor.processes import ProcessManager
from ucl.physiol.neuroconstruct.neuron import NeuronFileManager

# simulation control parameters
sim_path = '/home/eugenio/phd/nC_projects/if_network/'
project_filename = 'if_network.ncx'
sim_config_name = 'Default Simulation Configuration'
existing_stim = 'MF_stim'
nC_seed = 1234
simulator_seed = random.getrandbits(32)
active_mf_fraction = 0.3
sim_duration = 100.0
npatterns = 2
bias_values = [-0.005, -0.010]
# timestamp
now = time.localtime(time.time())
tstr = time.strftime("%Y-%m-%d.%H-%M-%S", now)

# load project and initialise
project_file = File(sim_path + project_filename)
print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))
pm = ProjectManager(None, None)
project = pm.loadProject(project_file)
print 'Loaded project: ' + project.getProjectName()

# generate network
pm.doGenerate(sim_config_name, nC_seed)
while pm.isGenerating():
    print 'Waiting for the project to be generated...'
    time.sleep(2)

# load existing simulation configurations and set sim duration and a couple of neuron options
sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(sim_duration)
project.neuronSettings.setNoConsole()
project.neuronSettings.setCopySimFiles(1)

# count generated cells
n_mf=project.generatedCellPositions.getNumberInCellGroup('MFs')
n_gr=project.generatedCellPositions.getNumberInCellGroup('GrCs')
n_generated_cells = project.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(n_generated_cells)

# simulate and record spiketimes
if n_generated_cells > 0:
    for n in range(npatterns):
        # delete all existing stimuli
        project.generatedElecInputs.reset()
        print "Preexisting stimuli: "+str(project.generatedElecInputs.getNumberSingleInputs(existing_stim))
        # generate random stimulation pattern (at fixed sparsity)
        active_inputs = random.sample(range(n_mf), int(round(n_mf*active_mf_fraction)))

        for bias in bias_values:
            # innermost loop: determine the simulation reference name
            sim_ref = tstr + "_p" + str(n) + "_b" + str(bias)
            project.simulationParameters.setReference(sim_ref)
            pattern_file=open(sim_path+"simulations/"+sim_ref+'_pattern.txt',"w")
            pattern_file.write(str(active_inputs))
            #Set the thresholding current
            bias_input = project.elecInputInfo.getStim("bias")
            bias_input.setAmp(NumberGenerator(bias))
            bias_input.setDur(NumberGenerator(sim_duration))
            project.elecInputInfo.updateStim(bias_input)
            # regenerate network structure, so that bias change can take effect
            pm.doGenerate(sim_config_name, nC_seed)
            while pm.isGenerating():
                print 'Waiting for the project to be regenerated...'
                time.sleep(1)
            
            for k, mf in enumerate(active_inputs):
                project.generatedElecInputs.addSingleInput(existing_stim,'RandomSpikeTrain','MFs',mf,0,0,None)
            
            # generate and compile neuron files
            print "Trial %d/%d. Generating NEURON scripts..." % (n+1, npatterns)
            project.neuronFileManager.generateTheNeuronFiles(sim_config, None, NeuronFileManager.RUN_HOC,simulator_seed)
            compile_process = ProcessManager(project.neuronFileManager.getMainHocFile())
            compile_success = compile_process.compileFileWithNeuron(0,0)
            # simulate
            if compile_success:
                print "...success."
                print "Simulating: simulation reference %s" % str(sim_ref)
                pm.doRunNeuron(sim_config)
                propsfile_path=sim_path+"simulations/"+sim_ref+"/simulation.props"
                while os.path.exists(propsfile_path)==0:
                    time.sleep(2)
                print("")
                pattern_file.close()
else:
    print "No cells generated."

System.exit(0)
