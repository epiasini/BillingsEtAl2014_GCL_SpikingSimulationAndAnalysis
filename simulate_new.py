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

sim_path = '/home/eugenio/phd/nC_projects/if_network/'
project_filename = 'if_network.ncx'
sim_config_name = 'Default Simulation Configuration'
existing_stim = 'MF_stim'
nC_seed = 1234
simulator_seed = random.randint(1, 9999)
active_mf_fraction = .3

sim_ref_base = 'first_test'

sim_duration = 100.0
npatterns = 2

# load project and initialise
project_file = File(sim_path + project_filename)
print ('Loading project from file: ' + project_file.getAbsolutePath() + ", exists: " + str(project_file.exists()))

pm = ProjectManager(None, None)
project = pm.loadProject(project_file)
print 'Loaded project: ' + project.getProjectName()

sim_config = project.simConfigInfo.getSimConfig(sim_config_name)
sim_config.setSimDuration(sim_duration)
pm.doGenerate(sim_config_name, nC_seed)
#sim_config = project.simConfigInfo.getSimConfig(sim_config_name) #needed just if generating NML file afterwards, I guess

while pm.isGenerating():
    print 'Waiting for the project to be generated...'
    time.sleep(2)

n_mf=project.generatedCellPositions.getNumberInCellGroup('MFs')
n_gr=project.generatedCellPositions.getNumberInCellGroup('GrCs')
n_generated_cells = project.generatedCellPositions.getNumberInAllCellGroups()
print "Number of cells generated: " + str(n_generated_cells)

project.neuronSettings.setNoConsole()
project.neuronSettings.setCopySimFiles(1)

now = time.localtime(time.time())
tstr = time.strftime("%Y-%m-%d.%H-%M-%S", now)

if n_generated_cells > 0:
    for n in range(npatterns):
        project.generatedElecInputs.reset()
        print "Preexisting stimuli: "+str(project.generatedElecInputs.getNumberSingleInputs(existing_stim))
        print "Trial %d/%d. Generating pattern and adding stimuli..." % (n+1, npatterns)
        sim_ref = tstr + "_" + str(n)
        project.simulationParameters.setReference(sim_ref)
        pattern_file=open(sim_path+"simulations/"+sim_ref+'_pattern.txt',"w")
        active_inputs = random.sample(range(n_mf), int(round(n_mf*active_mf_fraction)))
        pattern_file.write(str(active_inputs))
        for k, mf in enumerate(active_inputs):
            project.generatedElecInputs.addSingleInput(existing_stim,'RandomSpikeTrain','MFs',mf,0,0,None)
        print "Total stimuli: "+str(project.generatedElecInputs.getNumberSingleInputs(existing_stim))
        
        print "Trial %d/%d. Generating NEURON scripts..." % (n+1, npatterns)
        project.neuronFileManager.generateTheNeuronFiles(sim_config, None, NeuronFileManager.RUN_HOC,simulator_seed)
        compile_process = ProcessManager(project.neuronFileManager.getMainHocFile())
        compile_success = compile_process.compileFileWithNeuron(0,0)

        if compile_success:
            print "...success."
            print "Simulating: simulation reference %s" % str(sim_ref)
            pm.doRunNeuron(sim_config)
            propsfile_path=sim_path+"simulations/"+sim_ref+"/simulation.props"
            while os.path.exists(propsfile_path)==0:
                time.sleep(2)

        print("")
        pattern_file.close()

System.exit(0)
